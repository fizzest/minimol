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
  server = spawn(python, ["-m", "minimol", "--no-browser", "--port", "0"], { stdio: ["ignore", "pipe", "inherit"] });
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
    center: { x: box.x + box.width / 2, y: box.y + box.height / 2 },
    model: () => page.evaluate(() => window.__minimol.model()),
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
    async clickAtom(i) {
      const p = await h.atomXY(i);
      await page.mouse.click(p.x, p.y);
    },
    async drag(from, to) {
      await page.mouse.move(from.x, from.y);
      await page.mouse.down();
      await page.mouse.move((from.x + to.x) / 2, (from.y + to.y) / 2, { steps: 4 });
      await page.mouse.move(to.x, to.y, { steps: 4 });
      await page.mouse.up();
    },
    /** Wait until the molecule differs from `before` (by identity). */
    async waitForChange(before) {
      await page.waitForFunction((n) => window.__minimol.history.stack.length !== n || window.__minimol.history.pos !== n - 1, before);
    },
    async idle() {
      await page.waitForFunction(() => !window.__minimol.state.busy);
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

test("draw a molecule by clicking, edit it, add H, optimize, undo/redo, delete", async () => {
  const h = await openEditor();
  const { page } = h;

  // Click on empty space: a carbon appears.
  await page.mouse.click(h.center.x - 80, h.center.y);
  assert.equal((await h.model()).atoms.length, 1);

  // Click the carbon twice: grow a chain C-C-C.
  await h.clickAtom(0);
  await h.clickAtom(1);
  let m = await h.model();
  assert.equal(m.atoms.length, 3);
  assert.deepEqual(
    m.bonds.map((b) => [b.a, b.b, b.order]),
    [
      [0, 1, 1],
      [1, 2, 1],
    ],
  );
  assert.ok(Math.abs(dist(m, 0, 1) - 1.52) < 1e-9, "new bond uses covalent radii");

  // Choose oxygen and click the last carbon: it becomes O.
  await page.click('[data-element="O"]');
  await h.clickAtom(2);
  m = await h.model();
  assert.equal(m.atoms[2].el, "O");
  assert.equal(m.atoms.length, 3);

  // Click the C-C bond: single -> double -> triple -> single.
  const orders = [];
  for (let k = 0; k < 3; k++) {
    await page.mouse.click(...Object.values(await h.bondXY(0)));
    orders.push((await h.model()).bonds[0].order);
  }
  assert.deepEqual(orders, [2, 3, 1]);
  assert.equal((await h.model()).atoms.length, 3, "clicking a bond must not add atoms");

  // Add hydrogens: ethanol.
  await page.click("#btn-addh");
  await page.waitForFunction(() => window.__minimol.model().atoms.length === 9);
  await h.idle();
  assert.match(await h.text("#info"), /C2H6O · 9 atoms · charge 0/);
  await h.shot("01-ethanol-drawn");

  // Optimize with MMFF94.
  await page.selectOption("#method", "MMFF94");
  await page.click("#btn-opt");
  await page.waitForFunction(() => /kcal\/mol/.test(document.querySelector("#energy").textContent));
  m = await h.model();
  assert.ok(dist(m, 0, 1) > 1.49 && dist(m, 0, 1) < 1.55, `C-C ${dist(m, 0, 1)}`);
  assert.ok(dist(m, 1, 2) > 1.39 && dist(m, 1, 2) < 1.45, `C-O ${dist(m, 1, 2)}`);
  const mmffEnergy = await h.text("#energy");
  await h.shot("02-ethanol-mmff");

  // Optimize with g-xTB when it is installed.
  const methods = await page.$$eval("#method option:not([disabled])", (os) => os.map((o) => o.value));
  if (methods.includes("g-xTB")) {
    await page.selectOption("#method", "g-xTB");
    assert.equal(await page.isVisible("#uhf"), true, "spin field shown for g-xTB");
    await page.click("#btn-opt");
    await page.waitForFunction(() => / Eh/.test(document.querySelector("#energy").textContent), null, { timeout: 120000 });
    m = await h.model();
    assert.ok(dist(m, 1, 2) > 1.38 && dist(m, 1, 2) < 1.47, `g-xTB C-O ${dist(m, 1, 2)}`);
    await h.shot("03-ethanol-gxtb");
  } else {
    console.log("# g-xTB not installed; skipping that part");
  }
  const energyAfter = await h.text("#energy");

  // Undo: the energy no longer describes the structure and disappears.
  await page.keyboard.press("Control+z");
  assert.notEqual(await h.text("#energy"), energyAfter);
  // Redo brings back the exact optimized model, and with it its energy.
  await page.keyboard.press("Control+Shift+z");
  assert.equal(await h.text("#energy"), energyAfter);
  assert.notEqual(mmffEnergy, "");

  // Delete tool: removing the O also removes its hydrogen.
  await page.click('[data-tool="delete"]');
  await h.clickAtom(2);
  m = await h.model();
  assert.equal(m.atoms.length, 7);
  assert.match(await h.text("#info"), /C2H5 · 7 atoms/);
  assert.equal(await h.text("#energy"), "", "energy hidden after an edit");

  // Select tool + Delete key on a carbon removes it with its hydrogens.
  await page.click('[data-tool="select"]');
  await h.clickAtom(0);
  await page.keyboard.press("Delete");
  assert.match(await h.text("#info"), /CH2 · 3 atoms/);

  assert.deepEqual(h.errors, []);
  await page.close();
});

test("drag to bond two atoms and to grow a new atom; navigation never edits", async () => {
  const h = await openEditor();
  const { page } = h;
  await page.mouse.click(h.center.x - 120, h.center.y);
  await page.mouse.click(h.center.x + 120, h.center.y);
  let m = await h.model();
  assert.equal(m.atoms.length, 2);
  assert.equal(m.bonds.length, 0);

  // Drag atom 0 -> atom 1 makes a bond.
  await h.drag(await h.atomXY(0), await h.atomXY(1));
  m = await h.model();
  assert.deepEqual(m.bonds, [{ a: 0, b: 1, order: 1 }]);

  // Drag from atom 1 to empty space grows a new N at a standard bond length.
  await page.click('[data-element="N"]');
  const p1 = await h.atomXY(1);
  await h.drag(p1, { x: p1.x + 10, y: p1.y - 90 });
  m = await h.model();
  assert.equal(m.atoms.length, 3);
  assert.equal(m.atoms[2].el, "N");
  assert.deepEqual(m.bonds[1], { a: 1, b: 2, order: 1 });
  assert.ok(Math.abs(dist(m, 1, 2) - (0.76 + 0.71)) < 1e-9);

  // Rotating (drag on empty space), panning (right drag) and zooming change the view only.
  const before = JSON.stringify(m);
  const screenBefore = await h.atomXY(2);
  await h.drag({ x: h.center.x, y: h.center.y + 250 }, { x: h.center.x + 120, y: h.center.y + 200 });
  await page.mouse.move(h.center.x, h.center.y + 250);
  await page.mouse.down({ button: "right" });
  await page.mouse.move(h.center.x + 40, h.center.y + 230, { steps: 5 });
  await page.mouse.up({ button: "right" });
  await page.mouse.wheel(0, 300);
  assert.equal(JSON.stringify(await h.model()), before);
  const screenAfter = await h.atomXY(2);
  assert.ok(Math.hypot(screenAfter.x - screenBefore.x, screenAfter.y - screenBefore.y) > 5, "view changed");
  assert.equal((await page.evaluate(() => window.__minimol.history.stack.length)), 5, "only edits are in history");

  assert.deepEqual(h.errors, []);
  await page.close();
});

test("open, save, SMILES insert, errors and charges", async () => {
  const h = await openEditor();
  const { page } = h;
  const benzene = readFileSync(new URL("../data/benzene.xyz", import.meta.url));
  await page.setInputFiles("#file-input", { name: "benzene.xyz", mimeType: "text/plain", buffer: benzene });
  await page.waitForFunction(() => window.__minimol.model().atoms.length === 12);
  await h.idle();
  let m = await h.model();
  assert.equal(m.bonds.filter((b) => b.order === 2).length, 3);
  assert.match(await h.text("#message"), /Opened benzene\.xyz/);
  await h.shot("04-benzene");

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
  await page.waitForFunction(() => window.__minimol.model().atoms.length === 15);
  m = await h.model();
  let closest = Infinity;
  for (let i = 0; i < 12; i++) for (let j = 12; j < 15; j++) closest = Math.min(closest, dist(m, i, j));
  assert.ok(closest > 2, `inserted water overlaps benzene (${closest.toFixed(2)} A)`);

  // Charges: select the water O and press "-" (formal charge applies to the selection).
  await page.click('[data-tool="select"]');
  await page.keyboard.press("Escape");
  await h.clickAtom(12);
  await page.keyboard.press("-");
  assert.equal((await h.model()).atoms[12].charge, -1);
  assert.match(await h.text("#info"), /charge -1/);
  await h.shot("05-charge");

  // Measurements: two selected atoms show a distance.
  await h.clickAtom(0);
  await page.keyboard.down("Shift");
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
  await page.waitForFunction(() => window.__minimol.model().atoms.length === 14);
  await h.idle();

  // Hold the optimization response for a second.
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
