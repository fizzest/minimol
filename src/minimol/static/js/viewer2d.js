// 2D (skeletal, ChemDraw-style) view of a model with z = 0.
//
// Same interface as the 3D Viewer (setModel, pick, screenToModel, fit, pan,
// zoom, ...) so the editor can switch between them. Carbon atoms are bare
// vertices; other atoms get a label with their implicit hydrogens ("OH",
// "NH2"); charges are shown as superscripts; atoms over their valence get a
// red box, like ChemDraw.

import { color } from "./elements.js";

const LABEL_FONT = "600 15px system-ui, -apple-system, 'Segoe UI', sans-serif";
const SMALL_FONT = "600 10px system-ui, -apple-system, 'Segoe UI', sans-serif";
const MAX_SCALE = 60; // px per Angstrom when fitting small molecules
const LABEL_COLORS = { C: "#111", H: "#111", N: "#1d3fbf", O: "#c4231c", S: "#9a7d00", P: "#c25e00", F: "#3d8b1f", Cl: "#2f8f2f", Br: "#8a2323", I: "#6d1a7a" };

export class Viewer2D {
  constructor(container) {
    this.container = container;
    this.canvas = document.createElement("canvas");
    this.canvas.className = "view2d";
    container.appendChild(this.canvas);
    this.ctx = this.canvas.getContext("2d");
    this.model = { atoms: [], bonds: [] };
    this.selection = [];
    this.hover = null;
    this.ghost = null;
    this.scale = 40;
    this.offset = { x: 0, y: 0 }; // screen position (CSS px, relative to the canvas) of model (0, 0)
    this._frame = 0;
    this._resize();
    new ResizeObserver(() => this._resize()).observe(container);
  }

  // ------------------------------------------------------------ geometry

  _rect() {
    return this.canvas.getBoundingClientRect();
  }

  _toScreen(x, y) {
    return { x: this.offset.x + x * this.scale, y: this.offset.y - y * this.scale };
  }

  modelToClient(p) {
    const r = this._rect();
    const s = this._toScreen(p[0], p[1]);
    return { x: r.left + s.x, y: r.top + s.y };
  }

  screenToModel(clientX, clientY) {
    const r = this._rect();
    return [(clientX - r.left - this.offset.x) / this.scale, -(clientY - r.top - this.offset.y) / this.scale, 0];
  }

  centerPoint() {
    const r = this._rect();
    return this.screenToModel(r.left + r.width / 2, r.top + r.height / 2);
  }

  viewAxes() {
    return { right: [1, 0, 0], up: [0, 1, 0], toCamera: [0, 0, 1] };
  }

  // ---------------------------------------------------------- navigation

  fit(model) {
    // The container's size: a hidden canvas measures 0 x 0.
    const r = { width: this.container.clientWidth, height: this.container.clientHeight };
    if (!model.atoms.length) {
      this.scale = 40;
      this.offset = { x: r.width / 2, y: r.height / 2 };
      this.requestRender();
      return;
    }
    const xs = model.atoms.map((a) => a.x);
    const ys = model.atoms.map((a) => a.y);
    const [x0, x1, y0, y1] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
    const w = Math.max(x1 - x0, 1.5);
    const h = Math.max(y1 - y0, 1.5);
    this.scale = Math.min(MAX_SCALE, (r.width - 120) / w, (r.height - 120) / h);
    this.offset = { x: r.width / 2 - ((x0 + x1) / 2) * this.scale, y: r.height / 2 + ((y0 + y1) / 2) * this.scale };
    this.requestRender();
  }

  rotate() {
    /* a 2D drawing is not rotated with the mouse */
  }

  pan(dx, dy) {
    this.offset.x += dx;
    this.offset.y += dy;
    this.requestRender();
  }

  zoom(factor, clientX, clientY) {
    // Zoom about the pointer (or the centre): the point under it stays put.
    const r = this._rect();
    const px = clientX === undefined ? r.width / 2 : clientX - r.left;
    const py = clientY === undefined ? r.height / 2 : clientY - r.top;
    const s = Math.min(400, Math.max(4, this.scale / factor));
    this.offset.x = px - ((px - this.offset.x) * s) / this.scale;
    this.offset.y = py - ((py - this.offset.y) * s) / this.scale;
    this.scale = s;
    this.requestRender();
  }

  // ------------------------------------------------------------- picking

  _labelText(a, i) {
    if (a.el !== "C" || this._degree(i) === 0) return a.el;
    return null;
  }

  _degree(i) {
    let n = 0;
    for (const b of this.model.bonds) if (b.a === i || b.b === i) n++;
    return n;
  }

  pick(clientX, clientY) {
    const r = this._rect();
    const px = clientX - r.left;
    const py = clientY - r.top;
    let best = null;
    let bestD = Infinity;
    this.model.atoms.forEach((a, i) => {
      const s = this._toScreen(a.x, a.y);
      const d = Math.hypot(s.x - px, s.y - py);
      const reach = this._labelText(a, i) ? 11 : 8;
      if (d < reach && d < bestD) {
        best = { kind: "atom", index: i };
        bestD = d;
      }
    });
    if (best) return best;
    this.model.bonds.forEach((b, k) => {
      const p = this._toScreen(this.model.atoms[b.a].x, this.model.atoms[b.a].y);
      const q = this._toScreen(this.model.atoms[b.b].x, this.model.atoms[b.b].y);
      const vx = q.x - p.x;
      const vy = q.y - p.y;
      const len2 = vx * vx + vy * vy || 1;
      const t = Math.max(0, Math.min(1, ((px - p.x) * vx + (py - p.y) * vy) / len2));
      const d = Math.hypot(p.x + t * vx - px, p.y + t * vy - py);
      if (d < 6 && d < bestD) {
        best = { kind: "bond", index: k };
        bestD = d;
      }
    });
    return best;
  }

  atomsInRect(x0, y0, x1, y1) {
    const [l, rr] = [Math.min(x0, x1), Math.max(x0, x1)];
    const [t, bb] = [Math.min(y0, y1), Math.max(y0, y1)];
    const out = [];
    this.model.atoms.forEach((a, i) => {
      const p = this.modelToClient([a.x, a.y, 0]);
      if (p.x >= l && p.x <= rr && p.y >= t && p.y <= bb) out.push(i);
    });
    return out;
  }

  // ------------------------------------------------------------ state

  setModel(model) {
    this.model = model;
    this.requestRender();
  }

  setSelection(indices) {
    this.selection = indices.slice();
    this.requestRender();
  }

  setHover(pick) {
    this.hover = pick;
    this.requestRender();
  }

  setGhost(from, to) {
    this.ghost = from && to ? { from, to } : null;
    this.requestRender();
  }

  // ----------------------------------------------------------- drawing

  _resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(1, this.container.clientWidth);
    const h = Math.max(1, this.container.clientHeight);
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.canvas.style.width = `${w}px`;
    this.canvas.style.height = `${h}px`;
    this.dpr = dpr;
    this.requestRender();
  }

  requestRender() {
    if (this._frame) return;
    this._frame = requestAnimationFrame(() => {
      this._frame = 0;
      this._draw();
    });
  }

  /** Label pieces for an atom: [{text, sub}] plus where the H goes. */
  _label(a, i) {
    const sym = this._labelText(a, i);
    if (!sym) return null;
    const parts = [{ text: sym }];
    const h = a.h || 0;
    if (h > 0) {
      const hPart = [{ text: "H" }, ...(h > 1 ? [{ text: String(h), sub: true }] : [])];
      // Put the hydrogens on the side away from the bonds (ChemDraw writes "HO" for a bond on the right).
      let dx = 0;
      for (const b of this.model.bonds) {
        if (b.a === i) dx += this.model.atoms[b.b].x - a.x;
        else if (b.b === i) dx += this.model.atoms[b.a].x - a.x;
      }
      return dx > 0.3 ? { parts: [...hPart, ...parts], anchor: hPart.length } : { parts: [...parts, ...hPart], anchor: 0 };
    }
    return { parts, anchor: 0 };
  }

  _draw() {
    const g = this.ctx;
    const m = this.model;
    g.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    g.fillStyle = "#fff";
    g.fillRect(0, 0, this.canvas.width, this.canvas.height);
    const S = (i) => this._toScreen(m.atoms[i].x, m.atoms[i].y);
    const labels = m.atoms.map((a, i) => this._label(a, i));
    const line = Math.max(1.2, Math.min(2.2, this.scale * 0.04));

    // Highlights under the drawing.
    for (const i of this.selection) {
      if (!m.atoms[i]) continue;
      const p = S(i);
      g.fillStyle = "rgba(255, 176, 0, 0.45)";
      g.beginPath();
      g.arc(p.x, p.y, 11, 0, 2 * Math.PI);
      g.fill();
    }
    if (this.hover) {
      g.fillStyle = "rgba(30, 144, 255, 0.28)";
      g.strokeStyle = "rgba(30, 144, 255, 0.35)";
      if (this.hover.kind === "atom" && m.atoms[this.hover.index]) {
        const p = S(this.hover.index);
        g.beginPath();
        g.arc(p.x, p.y, 12, 0, 2 * Math.PI);
        g.fill();
      } else if (this.hover.kind === "bond" && m.bonds[this.hover.index]) {
        const b = m.bonds[this.hover.index];
        const p = S(b.a);
        const q = S(b.b);
        g.lineWidth = 12;
        g.lineCap = "round";
        g.beginPath();
        g.moveTo(p.x, p.y);
        g.lineTo(q.x, q.y);
        g.stroke();
      }
    }

    // Bonds.
    g.strokeStyle = "#111";
    g.lineWidth = line;
    g.lineCap = "round";
    m.bonds.forEach((b) => {
      let p = S(b.a);
      let q = S(b.b);
      const len = Math.hypot(q.x - p.x, q.y - p.y);
      if (len < 1e-6) return;
      const ux = (q.x - p.x) / len;
      const uy = (q.y - p.y) / len;
      // Leave room for labels.
      const gapA = labels[b.a] ? 9 : 0;
      const gapB = labels[b.b] ? 9 : 0;
      p = { x: p.x + ux * gapA, y: p.y + uy * gapA };
      q = { x: q.x - ux * gapB, y: q.y - uy * gapB };
      const nx = -uy;
      const ny = ux;
      const seg = (a, c, off, shrinkA = 0, shrinkB = 0) => {
        g.beginPath();
        g.moveTo(a.x + nx * off + ux * shrinkA, a.y + ny * off + uy * shrinkA);
        g.lineTo(c.x + nx * off - ux * shrinkB, c.y + ny * off - uy * shrinkB);
        g.stroke();
      };
      const sep = Math.max(3, this.scale * 0.11);
      if (b.order === 1) seg(p, q, 0);
      else if (b.order === 3) {
        seg(p, q, 0);
        seg(p, q, sep);
        seg(p, q, -sep);
      } else {
        // Double bond: draw the second line on the side of the other substituents
        // (inside a ring), shortened; centred if there is no preferred side.
        const side = this._doubleSide(b, nx, ny);
        if (side === 0) {
          seg(p, q, sep / 2);
          seg(p, q, -sep / 2);
        } else {
          seg(p, q, 0);
          const shrink = len * 0.12;
          seg(p, q, side * sep * 1.3, labels[b.a] ? 0 : shrink, labels[b.b] ? 0 : shrink);
        }
      }
    });

    // Labels.
    m.atoms.forEach((a, i) => {
      const lab = labels[i];
      const p = S(i);
      if (lab) this._drawLabel(lab, p, LABEL_COLORS[a.el] || color(a.el));
      if (a.charge) {
        const n = Math.abs(a.charge);
        const text = (n > 1 ? n : "") + (a.charge > 0 ? "+" : "−");
        g.font = SMALL_FONT;
        g.fillStyle = "#111";
        g.textAlign = "left";
        g.textBaseline = "middle";
        g.fillText(text, p.x + (lab ? 8 : 5), p.y - (lab ? 10 : 8));
      }
      if (a.bad) {
        g.strokeStyle = "#e11d1d";
        g.lineWidth = 1.5;
        g.strokeRect(p.x - 10, p.y - 10, 20, 20);
      }
    });

    if (this.ghost) {
      const p = this._toScreen(this.ghost.from[0], this.ghost.from[1]);
      const q = this._toScreen(this.ghost.to[0], this.ghost.to[1]);
      g.setLineDash([5, 4]);
      g.strokeStyle = "#777";
      g.lineWidth = 1.5;
      g.beginPath();
      g.moveTo(p.x, p.y);
      g.lineTo(q.x, q.y);
      g.stroke();
      g.setLineDash([]);
      g.beginPath();
      g.arc(q.x, q.y, 4, 0, 2 * Math.PI);
      g.stroke();
    }
  }

  _doubleSide(b, nx, ny) {
    const m = this.model;
    const pa = m.atoms[b.a];
    let sum = 0;
    for (const bd of m.bonds) {
      for (const [x, y] of [
        [b.a, b.b],
        [b.b, b.a],
      ]) {
        const other = bd.a === x && bd.b !== y ? bd.b : bd.b === x && bd.a !== y ? bd.a : -1;
        if (other < 0) continue;
        const o = m.atoms[other];
        // Screen y points down, model y up.
        sum += Math.sign((o.x - pa.x) * nx - (o.y - pa.y) * ny);
      }
    }
    return Math.sign(sum);
  }

  _drawLabel(lab, p, fill) {
    const g = this.ctx;
    // Measure each piece; the atom symbol (parts[anchor]) is centred on the atom.
    const widths = lab.parts.map((part) => {
      g.font = part.sub ? SMALL_FONT : LABEL_FONT;
      return g.measureText(part.text).width;
    });
    const before = widths.slice(0, lab.anchor).reduce((s, w) => s + w, 0);
    const total = widths.reduce((s, w) => s + w, 0);
    let x = p.x - before - widths[lab.anchor] / 2;
    g.fillStyle = "#fff";
    g.fillRect(x - 2, p.y - 9, total + 4, 18);
    g.fillStyle = fill;
    g.textAlign = "left";
    g.textBaseline = "middle";
    lab.parts.forEach((part, k) => {
      g.font = part.sub ? SMALL_FONT : LABEL_FONT;
      g.fillText(part.text, x, part.sub ? p.y + 5 : p.y);
      x += widths[k];
    });
  }
}
