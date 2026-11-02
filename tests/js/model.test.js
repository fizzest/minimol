import assert from "node:assert/strict";
import { test } from "node:test";

import { normalizeSymbol, SYMBOLS } from "../../src/minimol/static/js/elements.js";
import { History } from "../../src/minimol/static/js/history.js";
import * as M from "../../src/minimol/static/js/model.js";

const VIEW = { right: [1, 0, 0], up: [0, 1, 0], toCamera: [0, 0, 1] };

function deepFreeze(o) {
  Object.freeze(o);
  for (const v of Object.values(o)) if (v && typeof v === "object" && !Object.isFrozen(v)) deepFreeze(v);
  return o;
}

/** Ethanol heavy atoms C0-C1-O2 plus H3 on O2 and H4 on C0. */
function ethanolish() {
  return deepFreeze({
    atoms: [
      { el: "C", x: 0, y: 0, z: 0, charge: 0 },
      { el: "C", x: 1.5, y: 0, z: 0, charge: 0 },
      { el: "O", x: 2.0, y: 1.3, z: 0, charge: 0 },
      { el: "H", x: 2.9, y: 1.3, z: 0, charge: 0 },
      { el: "H", x: -0.5, y: -0.9, z: 0, charge: 0 },
    ],
    bonds: [
      { a: 0, b: 1, order: 1 },
      { a: 1, b: 2, order: 1 },
      { a: 2, b: 3, order: 1 },
      { a: 0, b: 4, order: 1 },
    ],
  });
}

test("element table has all 118 elements and normalises input", () => {
  assert.equal(SYMBOLS.length, 118);
  assert.equal(SYMBOLS[5], "C");
  assert.equal(SYMBOLS[117], "Og");
  assert.equal(normalizeSymbol(" cl "), "Cl");
  assert.equal(normalizeSymbol("FE"), "Fe");
  assert.equal(normalizeSymbol("Xx"), null);
  assert.equal(normalizeSymbol(""), null);
});

test("edits never mutate their input (inputs are deeply frozen)", () => {
  const m = ethanolish();
  M.addAtom(m, "N", [0, 0, 1]);
  M.addBond(m, 0, 2);
  M.setBondOrder(m, 0, 2);
  M.cycleBondOrder(m, 1);
  M.removeBond(m, 0);
  M.removeAtoms(m, [1]);
  M.removeAtomsWithHydrogens(m, [2]);
  M.setElement(m, [0], "N");
  M.changeCharge(m, [0], 1);
  M.append(m, m, [5, 0, 0]);
});

test("addAtom and addBond", () => {
  const { model, index } = M.addAtom(M.emptyModel(), "C", [1, 2, 3]);
  assert.equal(index, 0);
  assert.deepEqual(model.atoms[0], { el: "C", x: 1, y: 2, z: 3, charge: 0 });
  assert.throws(() => M.addAtom(model, "Qq", [0, 0, 0]));
  const m = ethanolish();
  assert.equal(M.addBond(m, 0, 1), m, "existing bond is a no-op");
  assert.equal(M.addBond(m, 1, 0), m, "existing bond in reverse is a no-op");
  assert.equal(M.addBond(m, 0, 0), m, "self bond is a no-op");
  assert.equal(M.addBond(m, 0, 99), m, "bad index is a no-op");
  const m2 = M.addBond(m, 0, 2, 2);
  assert.equal(m2.bonds.length, 5);
  assert.equal(M.validate(m2), null);
});

test("bond orders cycle 1 -> 2 -> 3 -> 1", () => {
  let m = ethanolish();
  const seen = [];
  for (let k = 0; k < 3; k++) {
    m = M.cycleBondOrder(m, 0);
    seen.push(m.bonds[0].order);
  }
  assert.deepEqual(seen, [2, 3, 1]);
  assert.equal(M.setBondOrder(m, 0, 1), m, "same order is a no-op");
});

test("removeAtoms reindexes bonds and returns an index map", () => {
  const { model, map } = M.removeAtoms(ethanolish(), [1]);
  assert.deepEqual(map, [0, -1, 1, 2, 3]);
  assert.equal(model.atoms.length, 4);
  assert.deepEqual(model.bonds, [
    { a: 1, b: 2, order: 1 },
    { a: 0, b: 3, order: 1 },
  ]);
  assert.equal(M.validate(model), null);
});

test("deleting an atom takes its dangling hydrogens with it, but not shared ones", () => {
  const m = ethanolish();
  assert.deepEqual(M.danglingHydrogens(m, [2]), [3]);
  const { model } = M.removeAtomsWithHydrogens(m, [2]);
  assert.equal(M.formula(model), "C2H");
  // An H bonded to two atoms (bridging) only goes if both go.
  const bridged = M.addBond(m, 3, 1);
  assert.deepEqual(M.danglingHydrogens(bridged, [2]), []);
  assert.deepEqual(M.danglingHydrogens(bridged, [1, 2]), [3]);
  // A lone H atom is not "dangling".
  const lone = M.addAtom(m, "H", [9, 9, 9]).model;
  assert.ok(!M.danglingHydrogens(lone, [0]).includes(5));
});

test("setElement and changeCharge return the same model when nothing changes", () => {
  const m = ethanolish();
  assert.equal(M.setElement(m, [0], "C"), m);
  assert.equal(M.setElement(m, [0], "N").atoms[0].el, "N");
  assert.throws(() => M.setElement(m, [0], "Zz"));
  let c = m;
  for (let k = 0; k < 20; k++) c = M.changeCharge(c, [0], +1);
  assert.equal(c.atoms[0].charge, M.MAX_CHARGE);
  assert.equal(M.changeCharge(c, [0], +1), c);
  assert.equal(M.totalCharge(M.changeCharge(m, [0, 2], -1)), -2);
});

test("formula uses Hill order", () => {
  assert.equal(M.formula(ethanolish()), "C2H2O");
  const salt = { atoms: [{ el: "Na" }, { el: "Cl" }], bonds: [] };
  assert.equal(M.formula(salt), "ClNa");
  assert.equal(M.formula(M.emptyModel()), "");
});

test("append shifts atoms and bond indices", () => {
  const m = ethanolish();
  const { model, firstIndex } = M.append(m, m, [10, 0, 0]);
  assert.equal(firstIndex, 5);
  assert.equal(model.atoms[5].x, 10);
  assert.deepEqual(model.bonds[4], { a: 5, b: 6, order: 1 });
  assert.equal(M.validate(model), null);
  const off = M.besideOffset(m, m);
  assert.ok(off[0] > 2.9 && Math.abs(off[1]) < 1e-12);
});

test("validate catches broken models", () => {
  const m = ethanolish();
  assert.equal(M.validate(m), null);
  assert.match(M.validate({ atoms: [{ el: "Qq", x: 0, y: 0, z: 0 }], bonds: [] }), /element/);
  assert.match(M.validate({ atoms: [{ el: "C", x: NaN, y: 0, z: 0 }], bonds: [] }), /coordinates/);
  assert.match(M.validate({ ...m, bonds: [...m.bonds, { a: 1, b: 0, order: 1 }] }), /duplicate/);
  assert.match(M.validate({ ...m, bonds: [{ a: 0, b: 9, order: 1 }] }), /bad atoms/);
  assert.match(M.validate({ ...m, bonds: [{ a: 0, b: 1, order: 1.5 }] }), /order/);
  assert.match(M.validate(null), /not a model/);
});

test("new atom directions point away from existing neighbours", () => {
  const single = { atoms: [{ el: "C", x: 0, y: 0, z: 0 }], bonds: [] };
  assert.deepEqual(M.newAtomDirection(single, 0, VIEW), [1, 0, 0]);

  const m = ethanolish();
  // Atom 1 has neighbours 0 and 2; the new direction must make > 90 deg with both.
  const d = M.newAtomDirection(m, 1, VIEW);
  for (const j of [0, 2]) {
    const u = M.vec.norm(M.vec.sub(M.pos(m.atoms[j]), M.pos(m.atoms[1])));
    assert.ok(M.vec.dot(d, u) < 0, `direction should point away from atom ${j}`);
  }
  assert.ok(Math.abs(M.vec.len(d) - 1) < 1e-9);

  // One neighbour: roughly tetrahedral angle to it.
  const one = M.newAtomDirection(m, 3, VIEW);
  const u = M.vec.norm(M.vec.sub(M.pos(m.atoms[2]), M.pos(m.atoms[3])));
  assert.ok(Math.abs((Math.acos(M.vec.dot(one, u)) * 180) / Math.PI - 109.47) < 0.01);

  // Linear neighbours (sum cancels): perpendicular to the axis.
  const linear = {
    atoms: [
      { el: "C", x: 0, y: 0, z: 0 },
      { el: "C", x: -1.5, y: 0, z: 0 },
      { el: "C", x: 1.5, y: 0, z: 0 },
    ],
    bonds: [
      { a: 0, b: 1, order: 1 },
      { a: 0, b: 2, order: 1 },
    ],
  };
  const p = M.newAtomDirection(linear, 0, VIEW);
  assert.ok(Math.abs(p[0]) < 1e-9 && Math.abs(M.vec.len(p) - 1) < 1e-9);

  // Trigonal planar (sum cancels): out of plane, towards the camera.
  const planar = {
    atoms: [
      { el: "B", x: 0, y: 0, z: 0 },
      { el: "H", x: 1, y: 0, z: 0 },
      { el: "H", x: -0.5, y: Math.sqrt(3) / 2, z: 0 },
      { el: "H", x: -0.5, y: -Math.sqrt(3) / 2, z: 0 },
    ],
    bonds: [1, 2, 3].map((b) => ({ a: 0, b, order: 1 })),
  };
  const z = M.newAtomDirection(planar, 0, VIEW);
  assert.ok(z[2] > 0.999);
});

test("new atom position uses covalent radii", () => {
  const single = { atoms: [{ el: "C", x: 0, y: 0, z: 0 }], bonds: [] };
  const p = M.newAtomPosition(single, 0, "H", VIEW);
  assert.ok(Math.abs(p[0] - (0.76 + 0.31)) < 1e-9);
});

test("measurements", () => {
  assert.equal(M.distance([0, 0, 0], [3, 4, 0]), 5);
  assert.ok(Math.abs(M.angle([1, 0, 0], [0, 0, 0], [0, 1, 0]) - 90) < 1e-9);
  assert.ok(Math.abs(M.angle([1, 0, 0], [0, 0, 0], [1, 0, 0])) < 1e-6);
  const dih = M.dihedral([1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]);
  assert.ok(Math.abs(Math.abs(dih) - 90) < 1e-9);
  assert.ok(Math.abs(M.dihedral([1, 0, 0], [0, 0, 0], [0, 0, 1], [1, 0, 1])) < 1e-9);
  assert.ok(Math.abs(Math.abs(M.dihedral([1, 0, 0], [0, 0, 0], [0, 0, 1], [-1, 0, 1])) - 180) < 1e-9);
  // Sign convention: opposite rotations give opposite signs.
  const a = M.dihedral([1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]);
  const b = M.dihedral([1, 0, 0], [0, 0, 0], [0, 0, 1], [0, -1, 1]);
  assert.ok(Math.abs(a + b) < 1e-9);
});

test("history: undo, redo, branch truncation and limit", () => {
  const h = new History("a", 3);
  assert.equal(h.push("a"), false);
  h.push("b");
  h.push("c");
  assert.equal(h.undo(), "b");
  assert.equal(h.undo(), "a");
  assert.equal(h.undo(), null);
  assert.equal(h.redo(), "b");
  h.push("x"); // drops "c"
  assert.equal(h.canRedo, false);
  assert.deepEqual(h.stack, ["a", "b", "x"]);
  h.push("y"); // limit 3 drops "a"
  assert.deepEqual(h.stack, ["b", "x", "y"]);
  assert.equal(h.current, "y");
  assert.equal(h.undo(), "x");
});
