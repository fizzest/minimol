// The molecule model: queries, validation and measurements.
//
// A model is { atoms: [{el, x, y, z, charge}], bonds: [{a, b, order}] }, the
// same JSON the server speaks. Models are never mutated: every edit (done by
// the server, see edit.py) produces a new model, so undo/redo is a plain list
// of models and the app can tell whether a result is stale by identity.

import { isElement } from "./elements.js";

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
