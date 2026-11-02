import assert from "node:assert/strict";
import { test } from "node:test";

import { normalizeSymbol, SYMBOLS } from "../../src/minimol/static/js/elements.js";
import { History } from "../../src/minimol/static/js/history.js";
import * as M from "../../src/minimol/static/js/model.js";

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

test("formula uses Hill order", () => {
  assert.equal(M.formula(ethanolish()), "C2H2O");
  const salt = { atoms: [{ el: "Na" }, { el: "Cl" }], bonds: [] };
  assert.equal(M.formula(salt), "ClNa");
  assert.equal(M.formula(M.emptyModel()), "");
});

test("model queries do not mutate their input", () => {
  const m = ethanolish(); // deeply frozen: any write would throw
  M.formula(m);
  M.neighbors(m, 1);
  M.findBond(m, 0, 1);
  M.totalCharge(m);
  M.append(m, m, [1, 0, 0]);
  M.besideOffset(m, m);
  assert.deepEqual(M.neighbors(m, 1), [0, 2]);
  assert.equal(M.findBond(m, 2, 1), 1);
  assert.equal(M.findBond(m, 0, 2), -1);
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
