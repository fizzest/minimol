// 3D view of a model: ball-and-stick rendering with thin outlines, picking,
// and trackball-style navigation.
//
// Scene graph: scene > pivot (rotation + pan) > content (shifted by -center)
// > meshes at plain model coordinates. The camera never rotates; it looks
// down -z, so the view plane is always the world xy-plane.

import * as THREE from "../vendor/three.module.js";
import { color, covalentRadius } from "./elements.js";

const BOND_RADIUS = 0.075;
const BOND_PICK_RADIUS = 0.17;
const OUTLINE = 0.028; // outline thickness in Angstrom
const FOV = 30;

export function atomRadius(el) {
  return el === "H" ? 0.2 : Math.min(0.55, 0.22 + 0.2 * covalentRadius(el));
}

const Y_AXIS = new THREE.Vector3(0, 1, 0);

export class Viewer {
  constructor(container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: false });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    container.appendChild(this.renderer.domElement);
    this.canvas = this.renderer.domElement;

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0xffffff);
    this.camera = new THREE.PerspectiveCamera(FOV, 1, 0.1, 2000);
    this.camera.position.set(0, 0, 20);

    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x8a8a8a, 1.4));
    const key = new THREE.DirectionalLight(0xffffff, 1.9);
    key.position.set(-4, 5, 10);
    this.scene.add(key);
    const fill = new THREE.DirectionalLight(0xffffff, 0.35);
    fill.position.set(5, -2, 4);
    this.scene.add(fill);

    this.pivot = new THREE.Group();
    this.content = new THREE.Group();
    this.pivot.add(this.content);
    this.scene.add(this.pivot);
    this.molGroup = new THREE.Group();
    this.overlayGroup = new THREE.Group();
    // Invisible, wider cylinders used only for picking bonds, so that a click
    // between the lines of a double bond still hits the bond.
    this.pickGroup = new THREE.Group();
    this.pickGroup.visible = false;
    this.content.add(this.molGroup, this.overlayGroup, this.pickGroup);

    this.sphereGeo = new THREE.SphereGeometry(1, 40, 28);
    this.cylGeo = new THREE.CylinderGeometry(1, 1, 1, 28, 1, false);
    this.outlineMat = new THREE.MeshBasicMaterial({ color: 0x1a1a1a, side: THREE.BackSide });
    this.materials = new Map();
    this.labelTextures = new Map();
    this.selectMat = new THREE.MeshBasicMaterial({ color: 0xffb000, transparent: true, opacity: 0.45, depthWrite: false });
    this.hoverMat = new THREE.MeshBasicMaterial({ color: 0x1e90ff, transparent: true, opacity: 0.35, depthWrite: false });
    this.ghostMat = new THREE.MeshBasicMaterial({ color: 0x555555, transparent: true, opacity: 0.55, depthWrite: false });
    this.pickMat = new THREE.MeshBasicMaterial();

    this.raycaster = new THREE.Raycaster();
    this.pickables = []; // what is drawn: atom spheres and bond cylinders
    this.bondPickers = []; // invisible, wider bond cylinders (fallback)
    this.model = { atoms: [], bonds: [] };
    this.selection = [];
    this.hover = null;
    this.ghost = null;
    this._frame = 0;

    this._resize();
    new ResizeObserver(() => this._resize()).observe(container);
  }

  // ------------------------------------------------------------ rendering

  requestRender() {
    if (this._frame) return;
    this._frame = requestAnimationFrame(() => {
      this._frame = 0;
      this.renderer.render(this.scene, this.camera);
    });
  }

  _resize() {
    const w = Math.max(1, this.container.clientWidth);
    const h = Math.max(1, this.container.clientHeight);
    this.renderer.setSize(w, h, true);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.requestRender();
  }

  _material(hex) {
    let m = this.materials.get(hex);
    if (!m) {
      m = new THREE.MeshPhongMaterial({ color: hex, specular: 0x2a2a2a, shininess: 40 });
      this.materials.set(hex, m);
    }
    return m;
  }

  _sphere(p, r, material, outline) {
    const mesh = new THREE.Mesh(this.sphereGeo, material);
    mesh.position.set(p[0], p[1], p[2]);
    mesh.scale.setScalar(r);
    if (outline) {
      const o = new THREE.Mesh(this.sphereGeo, this.outlineMat);
      o.position.copy(mesh.position);
      o.scale.setScalar(r + OUTLINE);
      this.molGroup.add(o);
    }
    return mesh;
  }

  _cylinder(p, q, r, material, outline) {
    const a = new THREE.Vector3(...p);
    const b = new THREE.Vector3(...q);
    const dir = b.clone().sub(a);
    const len = dir.length();
    if (len < 1e-6) return null;
    const quat = new THREE.Quaternion().setFromUnitVectors(Y_AXIS, dir.divideScalar(len));
    const mid = a.clone().add(b).multiplyScalar(0.5);
    const mesh = new THREE.Mesh(this.cylGeo, material);
    mesh.position.copy(mid);
    mesh.quaternion.copy(quat);
    mesh.scale.set(r, len, r);
    if (outline) {
      const o = new THREE.Mesh(this.cylGeo, this.outlineMat);
      o.position.copy(mid);
      o.quaternion.copy(quat);
      o.scale.set(r + OUTLINE, len, r + OUTLINE);
      this.molGroup.add(o);
    }
    return mesh;
  }

  _labelTexture(text) {
    let tex = this.labelTextures.get(text);
    if (!tex) {
      const c = document.createElement("canvas");
      c.width = c.height = 128;
      const g = c.getContext("2d");
      g.font = "bold 88px system-ui, sans-serif";
      g.textAlign = "center";
      g.textBaseline = "middle";
      g.lineWidth = 14;
      g.lineJoin = "round";
      g.strokeStyle = "rgba(255,255,255,0.9)";
      g.strokeText(text, 64, 64);
      g.fillStyle = "#111";
      g.fillText(text, 64, 64);
      tex = new THREE.CanvasTexture(c);
      tex.colorSpace = THREE.SRGBColorSpace;
      this.labelTextures.set(text, tex);
    }
    return tex;
  }

  /** Offsets of the parallel cylinders drawn for a bond of the given order. */
  _bondOffsets(m, bond) {
    const { a, b, order } = bond;
    if (order === 1) return { r: BOND_RADIUS, offsets: [0] };
    const pa = m.atoms[a];
    const pb = m.atoms[b];
    const axis = new THREE.Vector3(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z).normalize();
    // Put the bonds in the plane of a neighbouring atom when there is one.
    const others = [];
    for (const bd of m.bonds) {
      if (bd.a === a || bd.b === a || bd.a === b || bd.b === b) {
        for (const j of [bd.a, bd.b]) if (j !== a && j !== b) others.push(j);
      }
    }
    let ref = null;
    for (const j of others) {
      const o = m.atoms[j];
      const v = new THREE.Vector3(o.x - pa.x, o.y - pa.y, o.z - pa.z);
      v.sub(axis.clone().multiplyScalar(v.dot(axis)));
      if (v.lengthSq() > 1e-4) {
        ref = v.normalize();
        break;
      }
    }
    if (!ref) {
      ref = new THREE.Vector3().crossVectors(axis, Math.abs(axis.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0)).normalize();
    }
    const d = order === 2 ? 0.075 : 0.115;
    const offsets = order === 2 ? [-d, d] : [-d, 0, d];
    return { r: order === 2 ? 0.05 : 0.04, offsets, ref };
  }

  setModel(model) {
    this.model = model;
    this.molGroup.clear();
    this.pickGroup.clear();
    this.pickables = [];
    this.bondPickers = [];
    model.atoms.forEach((atom, i) => {
      const mesh = this._sphere([atom.x, atom.y, atom.z], atomRadius(atom.el), this._material(color(atom.el)), true);
      mesh.userData = { kind: "atom", index: i };
      this.molGroup.add(mesh);
      this.pickables.push(mesh);
      if (atom.charge) {
        const n = Math.abs(atom.charge);
        const text = (n > 1 ? n : "") + (atom.charge > 0 ? "+" : "−");
        const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: this._labelTexture(text), depthTest: false }));
        const r = atomRadius(atom.el);
        sprite.position.set(atom.x, atom.y, atom.z);
        // Put the glyph (drawn at the sprite's centre) just outside the atom's
        // upper-right edge: sprite.center is the anchor in units of the sprite size.
        const size = 0.6;
        const offset = 0.72 * r + 0.12;
        sprite.scale.setScalar(size);
        sprite.center.set(0.5 - offset / size, 0.5 - offset / size);
        sprite.renderOrder = 10;
        this.molGroup.add(sprite);
      }
    });
    model.bonds.forEach((bond, k) => {
      const pa = model.atoms[bond.a];
      const pb = model.atoms[bond.b];
      const { r, offsets, ref } = this._bondOffsets(model, bond);
      for (const off of offsets) {
        const shift = ref ? [ref.x * off, ref.y * off, ref.z * off] : [0, 0, 0];
        const p = [pa.x + shift[0], pa.y + shift[1], pa.z + shift[2]];
        const q = [pb.x + shift[0], pb.y + shift[1], pb.z + shift[2]];
        const mid = [(p[0] + q[0]) / 2, (p[1] + q[1]) / 2, (p[2] + q[2]) / 2];
        for (const [from, to, el] of [[p, mid, pa.el], [mid, q, pb.el]]) {
          const mesh = this._cylinder(from, to, r, this._material(color(el)), true);
          if (!mesh) continue;
          mesh.userData = { kind: "bond", index: k };
          this.molGroup.add(mesh);
          this.pickables.push(mesh);
        }
      }
      // Picker: the part of the bond outside both atom spheres.
      const ra = atomRadius(pa.el) * 0.9;
      const rb = atomRadius(pb.el) * 0.9;
      const len = Math.hypot(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
      if (len > ra + rb + 0.02) {
        const u = [(pb.x - pa.x) / len, (pb.y - pa.y) / len, (pb.z - pa.z) / len];
        const from = [pa.x + u[0] * ra, pa.y + u[1] * ra, pa.z + u[2] * ra];
        const to = [pb.x - u[0] * rb, pb.y - u[1] * rb, pb.z - u[2] * rb];
        const picker = this._cylinder(from, to, BOND_PICK_RADIUS, this.pickMat, false);
        picker.userData = { kind: "bond", index: k };
        this.pickGroup.add(picker);
        this.bondPickers.push(picker);
      }
    });
    this._rebuildOverlay();
  }

  setSelection(indices) {
    this.selection = indices.slice();
    this._rebuildOverlay();
  }

  setHover(pick) {
    const same = (pick && this.hover && pick.kind === this.hover.kind && pick.index === this.hover.index) || (!pick && !this.hover);
    if (same) return;
    this.hover = pick;
    this._rebuildOverlay();
  }

  /** Show a dashed "rubber band" from model point `from` to `to` while drawing; null hides it. */
  setGhost(from, to) {
    this.ghost = from && to ? { from, to } : null;
    this._rebuildOverlay();
  }

  _rebuildOverlay() {
    this.overlayGroup.clear();
    const m = this.model;
    for (const i of this.selection) {
      const a = m.atoms[i];
      if (!a) continue;
      const mesh = new THREE.Mesh(this.sphereGeo, this.selectMat);
      mesh.position.set(a.x, a.y, a.z);
      mesh.scale.setScalar(atomRadius(a.el) + 0.14);
      this.overlayGroup.add(mesh);
    }
    if (this.hover) {
      if (this.hover.kind === "atom" && m.atoms[this.hover.index]) {
        const a = m.atoms[this.hover.index];
        const mesh = new THREE.Mesh(this.sphereGeo, this.hoverMat);
        mesh.position.set(a.x, a.y, a.z);
        mesh.scale.setScalar(atomRadius(a.el) + 0.1);
        this.overlayGroup.add(mesh);
      } else if (this.hover.kind === "bond" && m.bonds[this.hover.index]) {
        const b = m.bonds[this.hover.index];
        const pa = m.atoms[b.a];
        const pb = m.atoms[b.b];
        const mesh = this._cylinder([pa.x, pa.y, pa.z], [pb.x, pb.y, pb.z], BOND_PICK_RADIUS, this.hoverMat, false);
        if (mesh) this.overlayGroup.add(mesh);
      }
    }
    if (this.ghost) {
      const mesh = this._cylinder(this.ghost.from, this.ghost.to, 0.06, this.ghostMat, false);
      if (mesh) this.overlayGroup.add(mesh);
      const end = new THREE.Mesh(this.sphereGeo, this.ghostMat);
      end.position.set(...this.ghost.to);
      end.scale.setScalar(0.22);
      this.overlayGroup.add(end);
    }
    this.requestRender();
  }

  // ------------------------------------------------------------ navigation

  /** Centre and zoom on the model; `keepRotation=false` also resets the orientation. */
  fit(model, keepRotation = false) {
    let center = [0, 0, 0];
    let radius = 0;
    if (model.atoms.length) {
      const n = model.atoms.length;
      center = model.atoms.reduce((s, a) => [s[0] + a.x / n, s[1] + a.y / n, s[2] + a.z / n], [0, 0, 0]);
      for (const a of model.atoms) {
        radius = Math.max(radius, Math.hypot(a.x - center[0], a.y - center[1], a.z - center[2]) + atomRadius(a.el));
      }
    }
    this.content.position.set(-center[0], -center[1], -center[2]);
    this.pivot.position.set(0, 0, 0);
    if (!keepRotation) this.pivot.quaternion.identity();
    const halfFov = THREE.MathUtils.degToRad(FOV / 2);
    const fitRadius = Math.max(radius + 0.8, 4);
    const aspectFactor = Math.min(1, this.camera.aspect);
    this.camera.position.z = Math.min(1500, fitRadius / Math.sin(halfFov) / aspectFactor);
    this.camera.far = this.camera.position.z + fitRadius * 4 + 100;
    this.camera.updateProjectionMatrix();
    this.requestRender();
  }

  rotate(dx, dy) {
    const angle = Math.hypot(dx, dy) * 0.008;
    if (angle === 0) return;
    const axis = new THREE.Vector3(dy, dx, 0).normalize();
    this.pivot.quaternion.premultiply(new THREE.Quaternion().setFromAxisAngle(axis, angle));
    this.requestRender();
  }

  _worldPerPixel() {
    const depth = this.camera.position.z - this.pivot.position.z;
    return (2 * depth * Math.tan(THREE.MathUtils.degToRad(FOV / 2))) / Math.max(1, this.container.clientHeight);
  }

  pan(dx, dy) {
    const s = this._worldPerPixel();
    this.pivot.position.x += dx * s;
    this.pivot.position.y -= dy * s;
    this.requestRender();
  }

  zoom(factor) {
    this.camera.position.z = THREE.MathUtils.clamp(this.camera.position.z * factor, 2, 1500);
    this.camera.far = Math.max(this.camera.far, this.camera.position.z + 100);
    this.camera.updateProjectionMatrix();
    this.requestRender();
  }

  // -------------------------------------------------- coordinates & picking

  _ndc(clientX, clientY) {
    const rect = this.canvas.getBoundingClientRect();
    return new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
  }

  /**
   * The atom or bond under a screen point: {kind: "atom"|"bond", index} or null.
   * Whatever is drawn under the cursor wins; only when the cursor is over
   * empty space do the wider invisible bond pickers count, so thin bonds and
   * the gap between the lines of a double bond are easy to hit.
   */
  pick(clientX, clientY) {
    this.pivot.updateMatrixWorld(true);
    this.raycaster.setFromCamera(this._ndc(clientX, clientY), this.camera);
    const hit =
      this.raycaster.intersectObjects(this.pickables, false)[0] ||
      this.raycaster.intersectObjects(this.bondPickers, false)[0];
    return hit ? { ...hit.object.userData } : null;
  }

  /** Model-space point under the cursor, on the view plane through model point `anchor`. */
  screenToModel(clientX, clientY, anchor) {
    this.pivot.updateMatrixWorld(true);
    const world = this.content.localToWorld(new THREE.Vector3(...anchor));
    const plane = new THREE.Plane(new THREE.Vector3(0, 0, 1), -world.z);
    this.raycaster.setFromCamera(this._ndc(clientX, clientY), this.camera);
    const hit = new THREE.Vector3();
    if (!this.raycaster.ray.intersectPlane(plane, hit)) return anchor.slice();
    const local = this.content.worldToLocal(hit);
    return [local.x, local.y, local.z];
  }

  /** Model point at the centre of rotation. */
  centerPoint() {
    this.pivot.updateMatrixWorld(true);
    const local = this.content.worldToLocal(this.pivot.getWorldPosition(new THREE.Vector3()));
    return [local.x, local.y, local.z];
  }

  /** Screen axes expressed in model coordinates. */
  viewAxes() {
    const inv = this.pivot.quaternion.clone().invert();
    const f = (x, y, z) => {
      const v = new THREE.Vector3(x, y, z).applyQuaternion(inv);
      return [v.x, v.y, v.z];
    };
    return { right: f(1, 0, 0), up: f(0, 1, 0), toCamera: f(0, 0, 1) };
  }

  /** Indices of atoms whose centres project inside a client-space rectangle. */
  atomsInRect(x0, y0, x1, y1) {
    const [left, right] = [Math.min(x0, x1), Math.max(x0, x1)];
    const [top, bottom] = [Math.min(y0, y1), Math.max(y0, y1)];
    const out = [];
    this.model.atoms.forEach((a, i) => {
      const p = this.modelToClient([a.x, a.y, a.z]);
      if (p.x >= left && p.x <= right && p.y >= top && p.y <= bottom) out.push(i);
    });
    return out;
  }

  /** Client (CSS pixel) coordinates of a model point. */
  modelToClient(p) {
    this.pivot.updateMatrixWorld(true);
    const v = this.content.localToWorld(new THREE.Vector3(...p)).project(this.camera);
    const rect = this.canvas.getBoundingClientRect();
    return { x: rect.left + ((v.x + 1) / 2) * rect.width, y: rect.top + ((1 - v.y) / 2) * rect.height };
  }
}
