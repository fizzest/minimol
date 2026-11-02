// The molecule model and every edit on it, as pure functions.
//
// A model is { atoms: [{el, x, y, z, charge}], bonds: [{a, b, order}] }, the
// same JSON the server speaks. Models are never mutated: each edit returns a
// new model (unchanged atoms/bonds are shared). That makes undo/redo a plain
// list of models and lets the app tell whether a result is stale by identity.

import { bondLength, isElement } from "./elements.js";

export const MAX_CHARGE = 8;

export function emptyModel() {
  return { atoms: [], bonds: [] };
}

// ---------------------------------------------------------------- vectors

export const vec = {
  sub: (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]],
  add: (p, q) => [p[0] + q[0], p[1] + q[1], p[2] + q[2]],
  scale: (p, s) => [p[0] * s, p[1] * s, p[2] * s],
  dot: (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2],
  cross: (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]],
  len: (p) => Math.hypot(p[0], p[1], p[2]),
  norm(p) {
    const l = Math.hypot(p[0], p[1], p[2]);
    return l > 1e-12 ? [p[0] / l, p[1] / l, p[2] / l] : [0, 0, 0];
  },
};

export function pos(atom) {
  return [atom.x, atom.y, atom.z];
}

// --------------------------------------------------------------- queries

export function findBond(m, a, b) {
  return m.bonds.findIndex((bd) => (bd.a === a && bd.b === b) || (bd.a === b && bd.b === a));
}

export function neighbors(m, i) {
  const out = [];
  for (const bd of m.bonds) {
    if (bd.a === i) out.push(bd.b);
    else if (bd.b === i) out.push(bd.a);
  }
  return out;
}

export function totalCharge(m) {
  return m.atoms.reduce((s, a) => s + (a.charge || 0), 0);
}

/** Molecular formula in Hill order (C, H, then alphabetical; alphabetical if no C). */
export function formula(m) {
  const counts = new Map();
  for (const a of m.atoms) counts.set(a.el, (counts.get(a.el) || 0) + 1);
  const keys = [...counts.keys()].sort();
  const order = counts.has("C")
    ? ["C", ...(counts.has("H") ? ["H"] : []), ...keys.filter((k) => k !== "C" && k !== "H")]
    : keys;
  return order.map((k) => k + (counts.get(k) > 1 ? counts.get(k) : "")).join("");
}

export function centroid(m) {
  if (m.atoms.length === 0) return [0, 0, 0];
  const s = m.atoms.reduce((acc, a) => [acc[0] + a.x, acc[1] + a.y, acc[2] + a.z], [0, 0, 0]);
  return vec.scale(s, 1 / m.atoms.length);
}

export function boundingRadius(m, center) {
  let r = 0;
  for (const a of m.atoms) r = Math.max(r, vec.len(vec.sub(pos(a), center)));
  return r;
}

/** Check a model's structure; returns an error string or null. */
export function validate(m) {
  if (!m || !Array.isArray(m.atoms) || !Array.isArray(m.bonds)) return "not a model";
  for (const [i, a] of m.atoms.entries()) {
    if (!a || !isElement(a.el)) return `atom ${i + 1}: bad element`;
    if (![a.x, a.y, a.z].every(Number.isFinite)) return `atom ${i + 1}: bad coordinates`;
    if (!Number.isInteger(a.charge ?? 0)) return `atom ${i + 1}: bad charge`;
  }
  const seen = new Set();
  for (const [k, b] of m.bonds.entries()) {
    const n = m.atoms.length;
    if (!b || !Number.isInteger(b.a) || !Number.isInteger(b.b)) return `bond ${k + 1}: bad atoms`;
    if (b.a < 0 || b.b < 0 || b.a >= n || b.b >= n || b.a === b.b) return `bond ${k + 1}: bad atoms`;
    if (![1, 2, 3].includes(b.order)) return `bond ${k + 1}: bad order`;
    const key = Math.min(b.a, b.b) * n + Math.max(b.a, b.b);
    if (seen.has(key)) return `bond ${k + 1}: duplicate`;
    seen.add(key);
  }
  return null;
}

// ------------------------------------------------------------------ edits

export function addAtom(m, el, p, charge = 0) {
  if (!isElement(el)) throw new Error(`Unknown element ${el}`);
  const atom = { el, x: p[0], y: p[1], z: p[2], charge };
  return { model: { atoms: [...m.atoms, atom], bonds: m.bonds }, index: m.atoms.length };
}

/** Add a bond; returns the same model if it already exists or is invalid. */
export function addBond(m, a, b, order = 1) {
  const n = m.atoms.length;
  if (a === b || a < 0 || b < 0 || a >= n || b >= n || findBond(m, a, b) >= 0) return m;
  return { atoms: m.atoms, bonds: [...m.bonds, { a, b, order }] };
}

export function setBondOrder(m, bi, order) {
  if (!m.bonds[bi] || m.bonds[bi].order === order) return m;
  const bonds = m.bonds.slice();
  bonds[bi] = { ...bonds[bi], order };
  return { atoms: m.atoms, bonds };
}

export function cycleBondOrder(m, bi) {
  if (!m.bonds[bi]) return m;
  return setBondOrder(m, bi, (m.bonds[bi].order % 3) + 1);
}

export function removeBond(m, bi) {
  if (!m.bonds[bi]) return m;
  return { atoms: m.atoms, bonds: m.bonds.filter((_, k) => k !== bi) };
}

/**
 * Hydrogens that would be left dangling if `indices` were deleted: H atoms
 * whose only neighbours are all in `indices`.
 */
export function danglingHydrogens(m, indices) {
  const doomed = new Set(indices);
  const extra = [];
  m.atoms.forEach((a, i) => {
    if (a.el !== "H" || doomed.has(i)) return;
    const nb = neighbors(m, i);
    if (nb.length > 0 && nb.every((j) => doomed.has(j))) extra.push(i);
  });
  return extra;
}

/**
 * Delete atoms (and their bonds). Returns { model, map } where map[old] is the
 * new index of each atom, or -1 if it was deleted.
 */
export function removeAtoms(m, indices) {
  const doomed = new Set(indices);
  if (doomed.size === 0) return { model: m, map: m.atoms.map((_, i) => i) };
  const map = new Array(m.atoms.length).fill(-1);
  const atoms = [];
  m.atoms.forEach((a, i) => {
    if (!doomed.has(i)) {
      map[i] = atoms.length;
      atoms.push(a);
    }
  });
  const bonds = [];
  for (const b of m.bonds) {
    if (map[b.a] >= 0 && map[b.b] >= 0) bonds.push({ a: map[b.a], b: map[b.b], order: b.order });
  }
  return { model: { atoms, bonds }, map };
}

/** Delete atoms together with hydrogens that only hang on them. */
export function removeAtomsWithHydrogens(m, indices) {
  return removeAtoms(m, [...indices, ...danglingHydrogens(m, indices)]);
}

export function setElement(m, indices, el) {
  if (!isElement(el)) throw new Error(`Unknown element ${el}`);
  let changed = false;
  const atoms = m.atoms.slice();
  for (const i of indices) {
    if (atoms[i] && atoms[i].el !== el) {
      atoms[i] = { ...atoms[i], el };
      changed = true;
    }
  }
  return changed ? { atoms, bonds: m.bonds } : m;
}

export function changeCharge(m, indices, delta) {
  let changed = false;
  const atoms = m.atoms.slice();
  for (const i of indices) {
    if (!atoms[i]) continue;
    const c = Math.max(-MAX_CHARGE, Math.min(MAX_CHARGE, (atoms[i].charge || 0) + delta));
    if (c !== atoms[i].charge) {
      atoms[i] = { ...atoms[i], charge: c };
      changed = true;
    }
  }
  return changed ? { atoms, bonds: m.bonds } : m;
}

/** Add all atoms and bonds of `other`, shifted by `offset`. */
export function append(m, other, offset = [0, 0, 0]) {
  const base = m.atoms.length;
  const atoms = m.atoms.concat(
    other.atoms.map((a) => ({ ...a, x: a.x + offset[0], y: a.y + offset[1], z: a.z + offset[2] })),
  );
  const bonds = m.bonds.concat(other.bonds.map((b) => ({ a: b.a + base, b: b.b + base, order: b.order })));
  return { model: { atoms, bonds }, firstIndex: base };
}

/** Offset that places `other` to the right (+x) of `m` with a gap, centred in y and z. */
export function besideOffset(m, other, gap = 2.5) {
  if (m.atoms.length === 0 || other.atoms.length === 0) return [0, 0, 0];
  const maxX = Math.max(...m.atoms.map((a) => a.x));
  const minX = Math.min(...other.atoms.map((a) => a.x));
  const c1 = centroid(m);
  const c2 = centroid(other);
  return [maxX - minX + gap, c1[1] - c2[1], c1[2] - c2[2]];
}

// ------------------------------------------------------ drawing geometry

const TETRAHEDRAL = (109.47 * Math.PI) / 180;

function perpendicular(u, hint) {
  // Component of `hint` perpendicular to unit vector `u`; any perpendicular if degenerate.
  let p = vec.sub(hint, vec.scale(u, vec.dot(hint, u)));
  if (vec.len(p) < 1e-3) {
    const other = Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    p = vec.cross(u, other);
  }
  return vec.norm(p);
}

/**
 * Direction for a new atom bonded to atom i, pointing away from its existing
 * neighbours. `view` gives the screen axes in model coordinates
 * ({right, up, toCamera}) so that new atoms appear where the user can see them.
 */
export function newAtomDirection(m, i, view) {
  const p = pos(m.atoms[i]);
  const dirs = neighbors(m, i)
    .map((j) => vec.norm(vec.sub(pos(m.atoms[j]), p)))
    .filter((d) => vec.len(d) > 0.5);
  if (dirs.length === 0) return vec.norm(view.right);
  if (dirs.length === 1) {
    const u = dirs[0];
    const w = perpendicular(u, view.up);
    return vec.norm(vec.add(vec.scale(u, Math.cos(TETRAHEDRAL)), vec.scale(w, Math.sin(TETRAHEDRAL))));
  }
  const sum = dirs.reduce((s, d) => vec.add(s, d), [0, 0, 0]);
  if (vec.len(sum) > 0.2) return vec.norm(vec.scale(sum, -1));
  // Neighbours cancel out (linear or planar-symmetric): go perpendicular to them.
  let n = vec.cross(dirs[0], dirs[1]);
  n = vec.len(n) > 0.1 ? vec.norm(n) : perpendicular(dirs[0], view.up);
  return vec.dot(n, view.toCamera) < 0 ? vec.scale(n, -1) : n;
}

/** Position for a new `el` atom bonded to atom i. */
export function newAtomPosition(m, i, el, view) {
  const d = newAtomDirection(m, i, view);
  return vec.add(pos(m.atoms[i]), vec.scale(d, bondLength(m.atoms[i].el, el)));
}

// ---------------------------------------------------------- measurements

export function distance(p, q) {
  return vec.len(vec.sub(p, q));
}

/** Angle p-q-r in degrees. */
export function angle(p, q, r) {
  const a = vec.norm(vec.sub(p, q));
  const b = vec.norm(vec.sub(r, q));
  return (Math.acos(Math.max(-1, Math.min(1, vec.dot(a, b)))) * 180) / Math.PI;
}

/** Dihedral p-q-r-s in degrees, in (-180, 180]. */
export function dihedral(p, q, r, s) {
  const b0 = vec.sub(q, p);
  const b1 = vec.sub(r, q);
  const b2 = vec.sub(s, r);
  const n1 = vec.cross(b0, b1);
  const n2 = vec.cross(b1, b2);
  const m1 = vec.cross(n1, vec.norm(b1));
  const x = vec.dot(n1, n2);
  const y = vec.dot(m1, n2);
  return (Math.atan2(y, x) * 180) / Math.PI;
}
