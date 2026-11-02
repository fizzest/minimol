// End-to-end tests: start the real server, drive the editor in Chromium with
// real mouse and keyboard events, and check the resulting molecule.
//
// Run with:  npm run e2e   (needs `playwright` with Chromium, and the Python
// package installed in .venv or the interpreter given by MINIMOL_PYTHON).

import assert from "node:assert/strict";
import { spawn } from "node:child_process";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { after, before, test } from "node:test";

const require = createRequire(import.meta.url);
const { chromium } = require("playwright");

const SHOTS = process.env.MINIMOL_SCREENSHOTS; // optional folder for screenshots
let server;
let url;
let browser;

before(async () => {
  const python = process.env.MINIMOL_PYTHON || ".venv/bin/python";
  server = spawn(python, ["-m", "minimol", "--browser", "--no-browser", "--port", "0"], { stdio: ["ignore", "pipe", "inherit"] });
  url = await new Promise((resolve, reject) => {
    let out = "";
    server.stdout.on("data", (d) => {
      out += d;
      const m = out.match(/running at (http:\/\/\S+)/);
      if (m) resolve(m[1]);
    });
    server.on("exit", (code) => reject(new Error(`server exited with ${code}`)));
    setTimeout(() => reject(new Error("server did not start")), 20000);
  });
  browser = await chromium.launch();
});

after(async () => {
  await browser?.close();
  server?.kill();
});

async function openEditor() {
  const page = await browser.newPage({ viewport: { width: 1280, height: 800 }, acceptDownloads: true });
  const errors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  // Chromium logs every HTTP error response; 400 is how the server reports
  // user errors (e.g. a bad SMILES), which the tests trigger on purpose.
  page.on("console", (m) => m.type() === "error" && !/status of 400/.test(m.text()) && errors.push(m.text()));
  await page.goto(url);
  await page.waitForFunction(() => document.querySelector("#method").options.length > 0);
  const box = await page.locator("#viewport canvas").boundingBox();
  const h = {
    page,
    errors,
    box,
    center: { x: box.x + box.width / 2, y: box.y + box.height / 2 },
    model: () => page.evaluate(() => window.__minimol.model()),
    formula: () => page.evaluate(() => window.__minimol.M.formula(window.__minimol.model())),
    atomXY: (i) =>
      page.evaluate((i) => {
        const w = window.__minimol;
        return w.viewer.modelToClient(w.M.pos(w.model().atoms[i]));
      }, i),
    async bondXY(k) {
      const m = await h.model();
      const b = m.bonds[k];
      const [p, q] = await Promise.all([h.atomXY(b.a), h.atomXY(b.b)]);
      return { x: (p.x + q.x) / 2, y: (p.y + q.y) / 2 };
    },
    /** Index of an atom whose centre is not covered by anything closer to the camera. */
    async visibleAtom(predicate) {
      const m = await h.model();
      for (let i = 0; i < m.atoms.length; i++) {
        if (!predicate(m.atoms[i], i, m)) continue;
        const p = await h.atomXY(i);
        const pick = await page.evaluate(({ x, y }) => window.__minimol.viewer.pick(x, y), p);
        if (pick && pick.kind === "atom" && pick.index === i) return i;
      }
      throw new Error("no visible atom matches");
    },
    async clickAtom(i) {
      const p = await h.atomXY(i);
      await page.mouse.click(p.x, p.y);
      await h.idle();
    },
    async hoverAtom(i) {
      const p = await h.atomXY(i);
      await page.mouse.move(p.x, p.y);
    },
    async drag(from, to, button = "left") {
      await page.mouse.move(from.x, from.y);
      await page.mouse.down({ button });
      await page.mouse.move((from.x + to.x) / 2, (from.y + to.y) / 2, { steps: 4 });
      await page.mouse.move(to.x, to.y, { steps: 4 });
      await page.mouse.up({ button });
      await h.idle();
    },
    async idle() {
      await page.waitForFunction(() => !window.__minimol.state.busy);
    },
    async waitFormula(f) {
      await page.waitForFunction((f) => window.__minimol.M.formula(window.__minimol.model()) === f, f);
      await h.idle();
    },
    text: (sel) => page.textContent(sel),
    async shot(name) {
      if (SHOTS) await page.screenshot({ path: `${SHOTS}/${name}.png` });
    },
  };
  return h;
}

function dist(m, i, j) {
  const a = m.atoms[i];
  const b = m.atoms[j];
  return Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
}

const carbons = (m) => m.atoms.map((a, i) => (a.el === "C" ? i : -1)).filter((i) => i >= 0);

test("drawing with hydrogens kept: grow, change element, bond orders, Ctrl+R, undo, right-click delete, box select", async () => {
  const h = await openEditor();
  const { page } = h;

  // Click on empty space: methane.
  await page.mouse.click(h.center.x - 80, h.center.y);
  await h.waitFormula("CH4");

  // Click the carbon: a hydrogen is replaced by CH3 -> ethane, then propane.
  await h.clickAtom(0);
  await h.waitFormula("C2H6");
  let m = await h.model();
  const c2 = carbons(m)[1];
  assert.ok(Math.abs(dist(m, 0, c2) - 1.52) < 1e-6);
  await h.clickAtom(c2);
  await h.waitFormula("C3H8");

  // Oxygen on the terminal carbon: that CH3 becomes OH -> ethanol.
  m = await h.model();
  const terminal = carbons(m)[2];
  await page.click('[data-element="O"]');
  await h.clickAtom(terminal);
  await h.waitFormula("C2H6O");

  // Click the C-C bond: double, triple, single, with hydrogens following.
  // (Re-deriving hydrogens renumbers atoms, so the bond is looked up each time.)
  const seen = [];
  for (let n = 0; n < 3; n++) {
    const mm = await h.model();
    const kk = mm.bonds.findIndex((b) => mm.atoms[b.a].el === "C" && mm.atoms[b.b].el === "C");
    const p = await h.bondXY(kk);
    await page.mouse.click(p.x, p.y);
    await h.idle();
    seen.push(await h.formula());
  }
  assert.deepEqual(seen, ["C2H4O", "C2H2O", "C2H6O"]);
  await h.shot("10-ethanol-drawn");

  // Ctrl+R optimizes (and must not reload the page).
  const historyLength = await page.evaluate(() => window.__minimol.history.stack.length);
  await page.selectOption("#method", "MMFF94");
  await page.mouse.move(h.center.x + 400, h.center.y + 250); // away from atoms
  await page.keyboard.press("Control+r");
  await page.waitForFunction(() => /kcal\/mol/.test(document.querySelector("#energy").textContent));
  assert.equal(await page.evaluate(() => window.__minimol.history.stack.length), historyLength + 1, "page was not reloaded");
  m = await h.model();
  const cc = m.bonds.filter((b) => m.atoms[b.a].el === "C" && m.atoms[b.b].el === "C");
  for (const b of cc) assert.ok(dist(m, b.a, b.b) > 1.49 && dist(m, b.a, b.b) < 1.55);

  // g-xTB, when installed.
  const methods = await page.$$eval("#method option:not([disabled])", (os) => os.map((o) => o.value));
  if (methods.includes("g-xTB")) {
    await page.selectOption("#method", "g-xTB");
    await page.click("#btn-opt");
    await page.waitForFunction(() => / Eh/.test(document.querySelector("#energy").textContent), null, { timeout: 120000 });
    await h.shot("11-ethanol-gxtb");
  }
  const energy = await h.text("#energy");
  await page.keyboard.press("Control+z");
  assert.notEqual(await h.text("#energy"), energy);
  await page.keyboard.press("Control+Shift+z");
  assert.equal(await h.text("#energy"), energy);

  // Right-click the oxygen in Draw mode: it is deleted and the carbon capped with H -> ethane.
  const o = await h.visibleAtom((a) => a.el === "O");
  const po = await h.atomXY(o);
  await page.mouse.click(po.x, po.y, { button: "right" });
  await h.waitFormula("C2H6");
  assert.equal(await h.text("#energy"), "", "energy hidden after an edit");

  // Ctrl+drag a box around everything selects every atom; Delete removes them.
  await page.mouse.move(h.box.x + 5, h.box.y + 5);
  await page.keyboard.down("Control");
  await h.drag({ x: h.box.x + 5, y: h.box.y + 5 }, { x: h.box.x + h.box.width - 5, y: h.box.y + h.box.height - 5 });
  await page.keyboard.up("Control");
  assert.equal((await page.evaluate(() => window.__minimol.state.selection.length)), 8);
  assert.match(await h.text("#measure"), /8 atoms selected/);
  await page.keyboard.press("Delete");
  await h.waitFormula("");

  assert.deepEqual(h.errors, []);
  await page.close();
});

test("fragment hotkeys: chain building, labels, nitro on a hydrogen, bond order on a bond", async () => {
  const h = await openEditor();
  const { page } = h;
  await page.fill("#smiles", "C");
  await page.press("#smiles", "Enter");
  await h.waitFormula("CH4");

  // Hover the carbon and type 1 1 o: each 1 adds a carbon and moves the target
  // to it; o turns that last carbon into OH -> ethanol.
  await h.hoverAtom(0);
  for (const key of ["1", "1", "o"]) {
    await page.keyboard.press(key);
    await h.idle();
  }
  assert.equal(await h.formula(), "C2H6O");
  let m = await h.model();
  const o = m.atoms.findIndex((a) => a.el === "O");
  const target = await page.evaluate(() => window.__minimol.state.target);
  assert.deepEqual(target, { kind: "atom", index: o }, "target moved to the new O");

  // A label key on a hydrogen turns that hydrogen into the group.
  const onCarbon = (i, mm) => mm.bonds.some((b) => (b.a === i && mm.atoms[b.b].el === "C") || (b.b === i && mm.atoms[b.a].el === "C"));
  const hIdx = await h.visibleAtom((a, i, mm) => a.el === "H" && onCarbon(i, mm));
  await h.hoverAtom(hIdx);
  await page.keyboard.press("Shift+N");
  await h.idle();
  assert.equal(await h.formula(), "C2H5NO3");
  m = await h.model();
  assert.equal(m.atoms.reduce((s, a) => s + a.charge, 0), 0, "nitro group is neutral overall");
  await h.shot("12-hotkeys");

  // Hover a C-C bond and press 2: double bond, hydrogens adjusted.
  const k = m.bonds.findIndex((b) => m.atoms[b.a].el === "C" && m.atoms[b.b].el === "C");
  const p = await h.bondXY(k);
  await page.mouse.move(p.x, p.y);
  await page.keyboard.press("2");
  await h.idle();
  m = await h.model();
  assert.equal(m.bonds.filter((b) => m.atoms[b.a].el === "C" && m.atoms[b.b].el === "C" && b.order === 2).length, 1);

  // With hotkeys switched off, keys over an atom do the plain thing (2 = Draw tool).
  await page.uncheck("#hotkeys");
  await page.click('[data-tool="select"]');
  await h.hoverAtom(0);
  const before = JSON.stringify(await h.model());
  await page.keyboard.press("2");
  assert.equal(JSON.stringify(await h.model()), before);
  assert.equal(await page.evaluate(() => window.__minimol.state.tool), "draw");

  assert.deepEqual(h.errors, []);
  await page.close();
});

test("skeleton drawing (Keep hydrogens off): drag to bond, drag to grow; navigation never edits", async () => {
  const h = await openEditor();
  const { page } = h;
  await page.uncheck("#keep-h");
  await page.mouse.click(h.center.x - 120, h.center.y);
  await h.idle();
  await page.mouse.click(h.center.x + 120, h.center.y);
  await h.idle();
  let m = await h.model();
  assert.equal(m.atoms.length, 2);
  assert.equal(m.bonds.length, 0);

  // Drag atom 0 -> atom 1 makes a bond.
  await h.drag(await h.atomXY(0), await h.atomXY(1));
  m = await h.model();
  assert.deepEqual(m.bonds, [{ a: 0, b: 1, order: 1 }]);

  // Drag from atom 1 to empty space grows an N along the drag, at a standard bond length.
  await page.click('[data-element="N"]');
  const p1 = await h.atomXY(1);
  await h.drag(p1, { x: p1.x + 10, y: p1.y - 90 });
  m = await h.model();
  assert.equal(m.atoms.length, 3);
  assert.equal(m.atoms[2].el, "N");
  assert.ok(Math.abs(dist(m, 1, 2) - (0.76 + 0.71)) < 1e-6);
  const shown = await h.atomXY(2);
  assert.ok(shown.y < p1.y - 30, "new atom is where it was dragged");

  // Rotate, pan and zoom change the view only.
  const before = JSON.stringify(m);
  const screenBefore = await h.atomXY(2);
  await h.drag({ x: h.center.x, y: h.center.y + 250 }, { x: h.center.x + 120, y: h.center.y + 200 });
  await h.drag({ x: h.center.x, y: h.center.y + 250 }, { x: h.center.x + 40, y: h.center.y + 230 }, "right");
  await page.mouse.wheel(0, 300);
  assert.equal(JSON.stringify(await h.model()), before);
  const screenAfter = await h.atomXY(2);
  assert.ok(Math.hypot(screenAfter.x - screenBefore.x, screenAfter.y - screenBefore.y) > 5, "view changed");
  assert.equal(await page.evaluate(() => window.__minimol.history.stack.length), 5, "only edits are in history");

  assert.deepEqual(h.errors, []);
  await page.close();
});

test("open, save, SMILES insert, errors and charges", async () => {
  const h = await openEditor();
  const { page } = h;
  const benzene = readFileSync(new URL("../data/benzene.xyz", import.meta.url));
  await page.setInputFiles("#file-input", { name: "benzene.xyz", mimeType: "text/plain", buffer: benzene });
  await h.waitFormula("C6H6");
  let m = await h.model();
  assert.equal(m.bonds.filter((b) => b.order === 2).length, 3);
  assert.match(await h.text("#message"), /Opened benzene\.xyz/);

  // Save as XYZ downloads a file with the same atoms.
  await page.selectOption("#save-format", "xyz");
  const [download] = await Promise.all([page.waitForEvent("download"), page.click("#btn-save")]);
  assert.equal(download.suggestedFilename(), "benzene.xyz");
  const saved = readFileSync(await download.path(), "utf8").split("\n");
  assert.equal(saved[0], "12");
  assert.equal(saved.filter((l) => l.startsWith("C ")).length, 6);

  // A bad SMILES shows an error and changes nothing.
  await page.fill("#smiles", "C1CC");
  await page.click("#smiles-form button");
  await page.waitForFunction(() => document.querySelector("#message").classList.contains("error"));
  assert.match(await h.text("#message"), /Could not parse SMILES/);
  assert.equal((await h.model()).atoms.length, 12);

  // A good SMILES is inserted next to the existing molecule without overlap.
  await page.fill("#smiles", "O");
  await page.press("#smiles", "Enter");
  await h.waitFormula("C6H8O");
  m = await h.model();
  let closest = Infinity;
  for (let i = 0; i < 12; i++) for (let j = 12; j < 15; j++) closest = Math.min(closest, dist(m, i, j));
  assert.ok(closest > 2, `inserted water overlaps benzene (${closest.toFixed(2)} A)`);

  // Charge: hover the water O and press "-" -> hydroxide (one H fewer).
  await h.hoverAtom(12);
  await page.keyboard.press("-");
  await h.waitFormula("C6H7O");
  m = await h.model();
  assert.equal(m.atoms[12].charge, -1);
  assert.match(await h.text("#info"), /charge -1/);

  // Measurements: two selected atoms show a distance.
  await page.click('[data-tool="select"]');
  await h.clickAtom(0);
  await page.keyboard.down("Shift"); // mouse.click has no modifiers option
  await h.clickAtom(1);
  await page.keyboard.up("Shift");
  assert.match(await h.text("#measure"), /^C1–C2 {2}1\.\d{4} Å$/);

  assert.deepEqual(h.errors, []);
  await page.close();
});

test("edits are refused while a calculation runs, so results never apply to a changed structure", async () => {
  const h = await openEditor();
  const { page } = h;
  await page.fill("#smiles", "CCCC");
  await page.press("#smiles", "Enter");
  await h.waitFormula("C4H10");

  // Hold the optimization response until released.
  let release;
  const held = new Promise((r) => (release = r));
  await page.route("**/api/optimize", async (route) => {
    await held;
    await route.continue();
  });
  await page.selectOption("#method", "MMFF94");
  await page.click("#btn-opt");
  await page.waitForFunction(() => window.__minimol.state.busy);

  const before = JSON.stringify(await h.model());
  await page.click('[data-tool="draw"]');
  await page.mouse.click(h.center.x + 300, h.center.y + 200); // would add an atom
  await page.keyboard.press("Control+z"); // would undo
  assert.equal(JSON.stringify(await h.model()), before);
  assert.match(await h.text("#message"), /Please wait/);
  assert.equal(await page.isDisabled("#btn-addh"), true);

  release();
  await page.waitForFunction(() => /kcal\/mol/.test(document.querySelector("#energy").textContent));
  assert.equal((await h.model()).atoms.length, 14);
  assert.deepEqual(h.errors, []);
  await page.close();
});
