// minimol editor: wires the toolbar, mouse tools and keyboard to the model,
// the viewers and the server.
//
// Rules that keep the editor predictable:
//  * history.current is the only molecule; every change goes through commit().
//    A 2D drawing is a model with `dim: 2`; whether the editor is in 2D or 3D
//    follows from the current model, so undo/redo also crosses mode switches.
//  * Structural edits are done by the server (RDKit knows valences), one at a
//    time. While a server call runs (state.busy) nothing may change the
//    molecule, so a result always applies to the structure it was made for.
//  * The energy shown belongs to one exact model object and disappears as
//    soon as the structure changes.

import { api } from "./api.js";
import { bondLength, normalizeSymbol } from "./elements.js";
import { History } from "./history.js";
import * as M from "./model.js";
import { Viewer } from "./viewer.js";
import { Viewer2D } from "./viewer2d.js";

const $ = (id) => document.getElementById(id);
const CLICK_SLOP = 4; // pixels a pointer may move and still count as a click
const FLAT_BOND = 1.5; // bond length of 2D drawings (layout.BOND_LENGTH on the server)
const QUICK_ELEMENTS = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"];
const KEY_ELEMENTS = { h: "H", c: "C", n: "N", o: "O", f: "F", p: "P", s: "S" };

const state = {
  tool: "draw",
  element: "C",
  selection: [], // ordered atom indices (order matters for angles/dihedrals)
  target: null, // atom/bond that hotkeys act on: under the mouse, or the last group added
  pointer: null, // last pointer position over the view
  busy: false,
  energy: null, // { model, text }
  fileName: "molecule",
  saved: null, // model last saved or opened
  keepH: loadSetting("keepH", true),
  hotkeys: loadSetting("hotkeys", true),
  fragments: new Map(), // key -> {name, mode}
  bondKeys: new Map(), // key -> {op, order | size, aromatic}
  lastSwitch: null, // { from, to }: switching back from `to` restores `from` exactly
};

function loadSetting(name, fallback) {
  try {
    const v = localStorage.getItem(`minimol.${name}`);
    return v === null ? fallback : v === "1";
  } catch {
    return fallback;
  }
}

function saveSetting(name, value) {
  try {
    localStorage.setItem(`minimol.${name}`, value ? "1" : "0");
  } catch {
    /* storage unavailable: setting lasts for this session only */
  }
}

const history = new History(M.emptyModel());
state.saved = history.current;

const viewers = {};
try {
  viewers[3] = new Viewer($("viewport"));
} catch (err) {
  window.minimolWebGLFailed = true;
  $("hint").textContent =
    `3D graphics (WebGL) could not start: ${err.message} ` +
    "In the desktop app, try starting it with: minimol --software-gl";
  throw err;
}
viewers[2] = new Viewer2D($("viewport"));
let viewer = viewers[3];
viewers[2].canvas.hidden = true;
for (const v of Object.values(viewers)) {
  v.canvas.tabIndex = -1; // focusable, so clicking the view takes keyboard focus from inputs
  v.setModel(history.current);
  v.fit(history.current);
}

const model = () => history.current;
const dimOf = (m) => (m.dim === 2 ? 2 : 3);
const dim = () => dimOf(model());
/** Tag a server molecule with the editor's dimension. */
const withDim = (m, d) => (d === 2 ? { ...m, dim: 2 } : m);

// ------------------------------------------------------------ messages

function setMessage(text, kind = "info") {
  const el = $("message");
  el.textContent = text || "";
  el.className = kind === "info" ? "" : kind;
  el.title = text || "";
}

function showNotes(notes, fallback) {
  if (notes && notes.length) setMessage(notes.join(" "), "warn");
  else setMessage(fallback);
}

// ------------------------------------------------------- model changes

function canEdit() {
  if (state.busy) {
    setMessage("Please wait for the current calculation to finish.", "warn");
    return false;
  }
  return true;
}

/** Make `next` the current molecule. `view`: null | "fit" | "reset". */
function commit(next, { view = null } = {}) {
  const problem = M.validate(next);
  if (problem) {
    setMessage(`Refusing an invalid structure (${problem}).`, "error");
    return false;
  }
  if (!history.push(next)) return false;
  refreshModel(view);
  return true;
}

function activateViewer() {
  const v = viewers[dim()];
  if (v === viewer) return false;
  viewer.setHover(null);
  viewer.canvas.hidden = true;
  viewer = v;
  viewer.canvas.hidden = false;
  return true;
}

function refreshModel(view = null) {
  const m = model();
  const switched = activateViewer();
  state.selection = state.selection.filter((i) => i < m.atoms.length);
  viewer.setModel(m);
  viewer.setSelection(state.selection);
  if (view === "reset" || switched) viewer.fit(m, false);
  else if (view === "fit") viewer.fit(m, true);
  setTarget(state.pointer && !drag ? viewer.pick(state.pointer.x, state.pointer.y) : null);
  updateUI();
}

function setTarget(pick) {
  const m = model();
  if (pick && ((pick.kind === "atom" && !m.atoms[pick.index]) || (pick.kind === "bond" && !m.bonds[pick.index]))) pick = null;
  state.target = pick;
  viewer.setHover(pick);
  $("hover").textContent = describe(pick);
}

function setSelection(indices) {
  state.selection = [...new Set(indices)].filter((i) => i >= 0 && i < model().atoms.length);
  viewer.setSelection(state.selection);
  updateUI();
}

function undo() {
  if (!canEdit() || !history.canUndo) return;
  history.undo();
  state.selection = [];
  refreshModel();
}

function redo() {
  if (!canEdit() || !history.canRedo) return;
  history.redo();
  state.selection = [];
  refreshModel();
}

// --------------------------------------------------- server round trips

/** Run a server call with editing locked. Quick calls show no spinner. */
async function busyCall(text, fn, { quiet = false } = {}) {
  if (state.busy) return null;
  state.busy = true;
  const show = () => {
    $("busy-text").textContent = text;
    $("busy").hidden = false;
    updateUI();
  };
  const timer = quiet ? setTimeout(show, 400) : (show(), null);
  try {
    return await fn();
  } catch (err) {
    setMessage(err.message, "error");
    return null;
  } finally {
    clearTimeout(timer);
    state.busy = false;
    $("busy").hidden = true;
    updateUI();
  }
}

/**
 * Apply one structural edit on the server (see edit.py). Returns the server's
 * result, or null if the edit failed (the message is shown).
 */
async function serverEdit(op) {
  if (!canEdit()) return null;
  const m = model();
  const d = dimOf(m);
  const body = { molecule: m, op, keep_h: state.keepH, view: viewer.viewAxes(), dim: d };
  const r = await busyCall("Editing…", () => api("/api/edit", body), { quiet: true });
  if (!r || model() !== m) return null;
  state.selection = state.selection.map((i) => r.map[i]).filter((i) => i >= 0);
  if (!commit(withDim(r.molecule, d))) return null;
  if (r.notes.length) setMessage(r.notes.join(" "), "warn");
  else setMessage("");
  return r;
}

/** Add a molecule from the server next to what is there (or as the whole scene). */
function insertMolecule(molecule, d) {
  const m = model();
  if (m.atoms.length === 0) return commit(withDim(molecule, d), { view: "reset" });
  return commit(withDim(M.append(m, molecule, M.besideOffset(m, molecule)).model, d), { view: "fit" });
}

async function openFile(file) {
  if (!file || !canEdit()) return;
  if (file.size > 50 * 1024 * 1024) {
    setMessage("That file is too large (over 50 MB).", "error");
    return;
  }
  const text = await file.text();
  const d = dim();
  const r = await busyCall(`Opening ${file.name}…`, () => api("/api/import", { name: file.name, text, dim: d }));
  if (!r) return;
  state.selection = [];
  if (commit(withDim(r.molecule, d), { view: "reset" })) {
    state.fileName = file.name.replace(/\.[^.]*$/, "") || "molecule";
    state.saved = model();
    updateUI();
    showNotes(r.notes, `Opened ${file.name}.`);
  }
}

let pendingSave = null;

async function save() {
  const m = model();
  if (m.atoms.length === 0) {
    setMessage("Nothing to save yet.", "warn");
    return;
  }
  const fmt = $("save-format").value;
  const r = await busyCall("Saving…", () =>
    api("/api/export", { molecule: m, format: fmt, title: state.fileName, dim: dimOf(m) }),
  );
  if (!r) return;
  const url = URL.createObjectURL(new Blob([r.text], { type: "text/plain" }));
  const a = document.createElement("a");
  a.href = url;
  a.download = `${state.fileName}.${fmt}`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
  if (window.minimolDesktop) {
    // The desktop window shows a save dialog and reports back through
    // minimolSaved / minimolSaveFailed below.
    pendingSave = m;
    return;
  }
  state.saved = m;
  updateUI();
  setMessage(`Saved ${a.download}.`);
}

// Called by the desktop window after its save dialog.
window.minimolSaved = (path) => {
  if (pendingSave) state.saved = pendingSave;
  pendingSave = null;
  const name = String(path).split(/[\\/]/).pop();
  state.fileName = name.replace(/\.[^.]*$/, "") || state.fileName;
  updateUI();
  setMessage(`Saved ${path}.`);
};
window.minimolSaveFailed = (reason) => {
  pendingSave = null;
  setMessage(reason ? `Not saved: ${reason}` : "Save cancelled.", reason ? "error" : "info");
};

async function insertSmiles(smiles) {
  smiles = smiles.trim();
  if (!smiles || !canEdit()) return;
  const d = dim();
  const r = await busyCall("Building structure…", () => api("/api/smiles", { smiles, dim: d }));
  if (!r) return;
  insertMolecule(r.molecule, d);
  $("smiles").value = "";
  $("smiles").blur(); // hand the keyboard back to the view
  showNotes(r.notes, `Inserted ${smiles}.`);
}

async function addHydrogens() {
  if (!canEdit() || dim() === 2) return;
  const m = model();
  if (m.atoms.length === 0) return;
  const r = await busyCall("Adding hydrogens…", () => api("/api/addh", { molecule: m }));
  if (!r || model() !== m) return;
  commit(r.molecule);
  setMessage(r.notes.join(" "));
}

async function optimize() {
  if (!canEdit()) return;
  const m = model();
  if (m.atoms.length === 0) {
    setMessage("Nothing to optimize yet.", "warn");
    return;
  }
  if (dimOf(m) === 2) {
    const r = await busyCall("Cleaning up…", () => api("/api/clean2d", { molecule: m }));
    if (r && model() === m && commit(withDim(r.molecule, 2))) setMessage("Cleaned up the drawing.");
    return;
  }
  const method = $("method").value;
  if (!method) {
    setMessage("No optimization method is available.", "error");
    return;
  }
  const body = { molecule: m, method, fixed: $("fix-selected").checked ? state.selection : [] };
  if (method === "g-xTB") {
    const u = $("uhf").value.trim();
    if (u !== "") {
      const n = Number(u);
      if (!Number.isInteger(n) || n < 0) {
        setMessage("Unpaired electrons must be a whole number ≥ 0 (or blank for automatic).", "error");
        return;
      }
      body.uhf = n;
    }
  }
  const t0 = performance.now();
  const r = await busyCall(`Optimizing with ${method}…`, () => api("/api/optimize", body));
  if (!r) return;
  if (model() !== m) {
    setMessage("The structure changed during the optimization; the result was discarded.", "warn");
    return;
  }
  if (!commit(r.molecule)) return;
  const eh = r.energy_unit === "Eh";
  const value = eh ? r.energy.toFixed(8) : r.energy.toFixed(4);
  state.energy = {
    model: model(),
    text: `${r.method}: E = ${value} ${eh ? "Eh" : "kcal/mol"}${r.converged ? "" : " (not converged)"}`,
  };
  updateUI();
  const secs = ((performance.now() - t0) / 1000).toFixed(1);
  showNotes(r.notes, `Optimized with ${r.method} in ${secs} s${r.converged ? "" : " (not converged)"}.`);
}

/** Switch between the 3D editor and the 2D drawing. */
async function switchMode(target) {
  if (!canEdit() || target === dim()) return;
  const m = model();
  // Switching straight back without editing restores the exact structure.
  if (state.lastSwitch && state.lastSwitch.to === m && dimOf(state.lastSwitch.from) === target) {
    const from = state.lastSwitch.from;
    state.lastSwitch = { from: m, to: from };
    state.selection = [];
    commit(from, { view: "reset" });
    setMessage(target === 2 ? "Back to the 2D drawing." : "Back to the 3D structure.");
    return;
  }
  if (m.atoms.length === 0) {
    commit(withDim(M.emptyModel(), target), { view: "reset" });
    return;
  }
  const r = await busyCall(target === 2 ? "Drawing in 2D…" : "Building 3D structure…", () =>
    api("/api/convert", { molecule: m, to: target }),
  );
  if (!r || model() !== m) return;
  const next = withDim(r.molecule, target);
  state.selection = [];
  if (commit(next, { view: "reset" })) {
    state.lastSwitch = { from: m, to: next };
    showNotes(r.notes, target === 2 ? "2D drawing: hydrogens are implicit." : "Built a 3D structure (ETKDG + force field).");
  }
}

// ------------------------------------------------------------ clipboard

async function copySelection({ cut = false } = {}) {
  const m = model();
  if (m.atoms.length === 0) {
    setMessage("Nothing to copy.", "warn");
    return;
  }
  const atoms = state.selection.length ? state.selection : null;
  let r;
  try {
    r = await api("/api/copy", { molecule: m, atoms, dim: dimOf(m) });
    if (r.target === "text") await navigator.clipboard.writeText(r.text);
  } catch (err) {
    setMessage(`Copy failed: ${err.message}`, "error");
    return;
  }
  setMessage(`Copied ${atoms ? `${atoms.length} atoms` : "the structure"} (${r.smiles}).`);
  if (cut && atoms) await serverEdit({ type: "delete", atoms });
}

async function paste() {
  if (!canEdit()) return;
  const d = dim();
  const body = { dim: d };
  if (!window.minimolDesktop) {
    try {
      body.text = await navigator.clipboard.readText();
    } catch (err) {
      setMessage(`The browser did not allow reading the clipboard (${err.message}).`, "error");
      return;
    }
  }
  const r = await busyCall("Pasting…", () => api("/api/paste", body));
  if (!r) return;
  insertMolecule(r.molecule, d);
  showNotes(r.notes, `Pasted from ${r.source}.`);
}

// ------------------------------------------------------------- edits

function deleteTarget(pick) {
  if (!pick) return;
  const m = model();
  if (pick.kind === "atom") serverEdit({ type: "delete", atoms: [pick.index] });
  else {
    const b = m.bonds[pick.index];
    serverEdit({ type: "delete", bonds: [[b.a, b.b]] });
  }
}

function cycleBond(k) {
  const b = model().bonds[k];
  serverEdit({ type: "bond_order", a: b.a, b: b.b, order: (b.order % 3) + 1 });
}

function changeCharge(atoms, delta) {
  if (atoms.length === 0) {
    setMessage("Hover an atom or select atoms first, then change the charge.", "warn");
    return;
  }
  serverEdit({ type: "charge", atoms, delta });
}

async function hotkey(key) {
  const t = state.target;
  if (!t) return;
  if (t.kind === "bond") {
    const spec = state.bondKeys.get(key);
    const b = model().bonds[t.index];
    if (spec.op === "order") await serverEdit({ type: "bond_order", a: b.a, b: b.b, order: spec.order });
    else await serverEdit({ type: "fuse", a: b.a, b: b.b, size: spec.size, aromatic: spec.aromatic });
    return;
  }
  if (key === "+" || key === "=") return changeCharge([t.index], +1);
  if (key === "-") return changeCharge([t.index], -1);
  const r = await serverEdit({ type: "fragment", atom: t.index, key });
  if (r && r.focus !== null) setTarget({ kind: "atom", index: r.focus });
}

// ------------------------------------------------------------- display

function atomLabel(i) {
  return `${model().atoms[i].el}${i + 1}`;
}

function describe(pick) {
  const m = model();
  if (!pick) return "";
  if (pick.kind === "atom") {
    const a = m.atoms[pick.index];
    const nb = M.neighbors(m, pick.index).length;
    return (
      `${a.el}${pick.index + 1} · ${nb} bond${nb === 1 ? "" : "s"}` +
      (a.h ? ` · ${a.h} H` : "") +
      (a.charge ? ` · charge ${a.charge > 0 ? "+" : ""}${a.charge}` : "") +
      (a.bad ? " · too many bonds" : "")
    );
  }
  const b = m.bonds[pick.index];
  const length = dimOf(m) === 2 ? "" : ` · ${M.distance(M.pos(m.atoms[b.a]), M.pos(m.atoms[b.b])).toFixed(3)} Å`;
  return `bond ${atomLabel(b.a)}–${atomLabel(b.b)} · order ${b.order}${length}`;
}

function measurementText() {
  const m = model();
  const s = state.selection;
  if (dimOf(m) === 2) return s.length ? `${s.length} atom${s.length === 1 ? "" : "s"} selected` : "";
  const p = (i) => M.pos(m.atoms[i]);
  switch (s.length) {
    case 0:
      return "";
    case 1: {
      const [x, y, z] = p(s[0]);
      return `${atomLabel(s[0])} (${x.toFixed(3)}, ${y.toFixed(3)}, ${z.toFixed(3)})`;
    }
    case 2:
      return `${s.map(atomLabel).join("–")}  ${M.distance(p(s[0]), p(s[1])).toFixed(4)} Å`;
    case 3:
      return `∠ ${s.map(atomLabel).join("–")}  ${M.angle(p(s[0]), p(s[1]), p(s[2])).toFixed(2)}°`;
    case 4:
      return `dihedral ${s.map(atomLabel).join("–")}  ${M.dihedral(...s.map(p)).toFixed(2)}°`;
    default:
      return `${s.length} atoms selected`;
  }
}

function hasUnsavedChanges() {
  return model() !== state.saved && model().atoms.length > 0;
}

function updateUI() {
  const m = model();
  const flat = dimOf(m) === 2;
  const n = m.atoms.length;
  const q = M.totalCharge(m);
  $("info").textContent = n ? `${M.formula(m)} · ${n} atom${n === 1 ? "" : "s"} · charge ${q > 0 ? "+" : ""}${q}` : "Empty";
  $("measure").textContent = measurementText();
  $("energy").textContent = state.energy && state.energy.model === m ? state.energy.text : "";
  $("hint").hidden = n > 0;
  document.title = `${state.fileName}${hasUnsavedChanges() ? " •" : ""} — minimol`;

  for (const b of document.querySelectorAll("[data-tool]")) b.classList.toggle("active", b.dataset.tool === state.tool);
  for (const b of document.querySelectorAll("[data-element]")) b.classList.toggle("active", b.dataset.element === state.element);
  $("element-other").classList.toggle("active", !QUICK_ELEMENTS.includes(state.element));
  $("viewport").className = `tool-${state.tool}`;
  $("keep-h").checked = state.keepH;
  $("hotkeys").checked = state.hotkeys;
  $("mode-3d").classList.toggle("active", !flat);
  $("mode-2d").classList.toggle("active", flat);

  const locked = state.busy && !$("busy").hidden;
  for (const el of document.querySelectorAll("[data-edit], #method, #uhf, #fix-selected, #smiles, #element-other, #keep-h")) {
    el.disabled = locked;
  }
  if (!locked) {
    $("btn-undo").disabled = !history.canUndo;
    $("btn-redo").disabled = !history.canRedo;
  }
  // 2D drawings: hydrogens are implicit and "optimize" means "clean up".
  $("btn-addh").disabled = locked || flat;
  $("keep-h").disabled = locked || flat;
  $("method").disabled = locked || flat;
  $("fix-selected").disabled = locked || flat;
  $("btn-opt").textContent = flat ? "Clean up" : "Optimize";
  $("btn-opt").title = flat ? "Lay out the drawing again (Ctrl+R)" : "Optimize the geometry (Ctrl+R)";
  for (const opt of $("save-format").options) opt.disabled = flat && (opt.value === "xyz" || opt.value === "pdb");
  if ($("save-format").selectedOptions[0]?.disabled) $("save-format").value = "mol";
  $("uhf-wrap").hidden = flat || $("method").value !== "g-xTB";
}

function setTool(tool) {
  state.tool = tool;
  updateUI();
}

function setElement(el) {
  state.element = el;
  state.tool = "draw";
  if (QUICK_ELEMENTS.includes(el)) $("element-other").value = "";
  updateUI();
}

// ------------------------------------------------------- mouse tools

/** Where a new atom goes when dragging from atom i to a screen point. */
function dragTarget(i, clientX, clientY) {
  const a = model().atoms[i];
  const p = M.pos(a);
  const d = M.vec.sub(viewer.screenToModel(clientX, clientY, p), p);
  if (M.vec.len(d) < 1e-6) return null;
  if (dim() === 2) {
    // ChemDraw-style: bonds snap to 30 degree steps and have one length.
    const angle = Math.round(Math.atan2(d[1], d[0]) / (Math.PI / 6)) * (Math.PI / 6);
    return [p[0] + FLAT_BOND * Math.cos(angle), p[1] + FLAT_BOND * Math.sin(angle), 0];
  }
  return M.vec.add(p, M.vec.scale(M.vec.norm(d), bondLength(a.el, state.element)));
}

function handleClick(pick, e) {
  const m = model();
  if (state.tool === "select") {
    const additive = e.shiftKey;
    if (!pick) {
      if (!additive) setSelection([]);
    } else if (pick.kind === "atom") {
      const i = pick.index;
      if (additive) toggleSelected([i]);
      else setSelection([i]);
    } else {
      const b = m.bonds[pick.index];
      setSelection(additive ? [...state.selection, b.a, b.b] : [b.a, b.b]);
    }
    return;
  }
  if (state.tool === "draw") {
    if (!pick) {
      serverEdit({ type: "add_atom", element: state.element, pos: viewer.screenToModel(e.clientX, e.clientY, viewer.centerPoint()) });
    } else if (pick.kind === "atom") {
      const i = pick.index;
      if (m.atoms[i].el !== state.element) serverEdit({ type: "element", atom: i, element: state.element });
      else serverEdit({ type: "grow", atom: i, element: state.element });
    } else {
      cycleBond(pick.index);
    }
  } else if (state.tool === "delete") {
    deleteTarget(pick);
  }
}

function toggleSelected(indices) {
  let sel = state.selection.slice();
  for (const i of indices) sel = sel.includes(i) ? sel.filter((j) => j !== i) : [...sel, i];
  setSelection(sel);
}

function handleBondDrag(i, e) {
  const m = model();
  const target = viewer.pick(e.clientX, e.clientY);
  if (target && target.kind === "atom") {
    const j = target.index;
    if (j === i) return;
    const k = M.findBond(m, i, j);
    if (k >= 0) cycleBond(k);
    else serverEdit({ type: "bond", a: i, b: j });
    return;
  }
  const p = dragTarget(i, e.clientX, e.clientY);
  if (!p) return;
  serverEdit({ type: "grow", atom: i, element: state.element, hint: M.vec.norm(M.vec.sub(p, M.pos(m.atoms[i]))) });
}

let drag = null;

function onPointerDown(e) {
  const canvas = e.currentTarget;
  if (drag) return;
  canvas.focus({ preventScroll: true });
  canvas.setPointerCapture(e.pointerId);
  const pick = viewer.pick(e.clientX, e.clientY);
  let mode = dim() === 2 ? "box" : "rotate"; // in 2D, dragging on empty space selects
  if (e.button === 0 && (e.ctrlKey || e.metaKey)) mode = "box";
  else if (e.button === 1 || e.button === 2 || e.shiftKey) mode = "pan";
  else if (state.tool === "draw" && pick && pick.kind === "atom" && !state.busy) mode = "bond";
  drag = {
    id: e.pointerId,
    canvas,
    button: e.button,
    x0: e.clientX,
    y0: e.clientY,
    x: e.clientX,
    y: e.clientY,
    mode,
    pick,
    moved: false,
    ctrl: e.ctrlKey || e.metaKey,
  };
}

function onPointerMove(e) {
  state.pointer = { x: e.clientX, y: e.clientY };
  if (!drag) {
    setTarget(viewer.pick(e.clientX, e.clientY));
    return;
  }
  if (e.pointerId !== drag.id) return;
  if (!drag.moved && Math.hypot(e.clientX - drag.x0, e.clientY - drag.y0) > CLICK_SLOP) drag.moved = true;
  const dx = e.clientX - drag.x;
  const dy = e.clientY - drag.y;
  drag.x = e.clientX;
  drag.y = e.clientY;
  if (!drag.moved) return;
  if (drag.mode === "rotate") viewer.rotate(dx, dy);
  else if (drag.mode === "pan") viewer.pan(dx, dy);
  else if (drag.mode === "box") showBox(drag.x0, drag.y0, e.clientX, e.clientY);
  else if (drag.mode === "bond") {
    const i = drag.pick.index;
    const from = M.pos(model().atoms[i]);
    const target = viewer.pick(e.clientX, e.clientY);
    const onAtom = target && target.kind === "atom" && target.index !== i;
    viewer.setHover(onAtom ? target : null);
    const to = onAtom ? M.pos(model().atoms[target.index]) : dragTarget(i, e.clientX, e.clientY);
    viewer.setGhost(from, to);
  }
}

function showBox(x0, y0, x1, y1) {
  const rect = $("viewport").getBoundingClientRect();
  const box = $("selbox");
  box.hidden = false;
  box.style.left = `${Math.min(x0, x1) - rect.left}px`;
  box.style.top = `${Math.min(y0, y1) - rect.top}px`;
  box.style.width = `${Math.abs(x1 - x0)}px`;
  box.style.height = `${Math.abs(y1 - y0)}px`;
}

function endDrag(e, cancelled) {
  if (!drag || e.pointerId !== drag.id) return;
  const d = drag;
  drag = null;
  if (d.canvas.hasPointerCapture(e.pointerId)) d.canvas.releasePointerCapture(e.pointerId);
  viewer.setGhost(null);
  $("selbox").hidden = true;
  if (cancelled) return;
  if (d.mode === "box" && d.moved) {
    const inside = viewer.atomsInRect(d.x0, d.y0, e.clientX, e.clientY);
    setSelection(e.shiftKey ? [...state.selection, ...inside] : inside);
    return;
  }
  if (!d.moved && d.button === 0 && d.ctrl) {
    // Ctrl+click toggles an atom (or both atoms of a bond) in the selection.
    if (d.pick) {
      const m = model();
      toggleSelected(d.pick.kind === "atom" ? [d.pick.index] : [m.bonds[d.pick.index].a, m.bonds[d.pick.index].b]);
    }
    return;
  }
  if (!d.moved && d.button === 0) handleClick(d.pick, e);
  else if (!d.moved && d.button === 2 && state.tool !== "select") deleteTarget(d.pick);
  else if (d.moved && d.mode === "bond") handleBondDrag(d.pick.index, e);
}

for (const v of Object.values(viewers)) {
  const canvas = v.canvas;
  canvas.addEventListener("pointerdown", onPointerDown);
  canvas.addEventListener("pointermove", onPointerMove);
  canvas.addEventListener("pointerup", (e) => endDrag(e, false));
  canvas.addEventListener("pointercancel", (e) => endDrag(e, true));
  canvas.addEventListener("pointerleave", () => {
    state.pointer = null;
    if (!drag) setTarget(null);
  });
  canvas.addEventListener("contextmenu", (e) => e.preventDefault());
  canvas.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      const scale = e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? 400 : 1;
      viewer.zoom(Math.exp(Math.max(-200, Math.min(200, e.deltaY * scale)) * 0.0015), e.clientX, e.clientY);
    },
    { passive: false },
  );
}

// Drag and drop a file onto the view.
const viewport = $("viewport");
viewport.addEventListener("dragover", (e) => {
  e.preventDefault();
  viewport.classList.add("drop");
});
viewport.addEventListener("dragleave", () => viewport.classList.remove("drop"));
viewport.addEventListener("drop", (e) => {
  e.preventDefault();
  viewport.classList.remove("drop");
  openFile(e.dataTransfer.files[0]);
});

// ---------------------------------------------------------- keyboard

function hotkeyApplies(key) {
  if (!state.hotkeys || !state.target) return false;
  if (state.target.kind === "bond") return state.bondKeys.has(key);
  return state.fragments.has(key) || key === "+" || key === "=" || key === "-";
}

window.addEventListener("keydown", (e) => {
  const t = e.target;
  const typing = t instanceof HTMLElement && (t.tagName === "INPUT" || t.tagName === "TEXTAREA") && t.type !== "checkbox";
  const mod = e.ctrlKey || e.metaKey;
  const key = e.key.toLowerCase();
  if (typing) {
    if (e.key === "Escape") t.blur();
    return;
  }
  if (t instanceof HTMLElement && t.tagName === "SELECT" && !mod) {
    // Arrow keys etc. belong to the drop-down; give other keys back to the editor.
    if (e.key.startsWith("Arrow") || e.key === "Enter" || e.key === " ") return;
    t.blur();
  }
  if (mod) {
    const actions = {
      z: () => (e.shiftKey ? redo() : undo()),
      y: redo,
      s: save,
      r: optimize,
      c: () => copySelection(),
      x: () => copySelection({ cut: true }),
      v: paste,
      o: () => canEdit() && $("file-input").click(),
      a: () => setSelection(model().atoms.map((_, i) => i)),
    };
    if (actions[key] && !e.altKey) {
      e.preventDefault();
      actions[key]();
    }
    return;
  }
  if (e.altKey) return;
  if (e.key === "Tab") {
    e.preventDefault();
    switchMode(dim() === 2 ? 3 : 2);
    return;
  }
  if (e.key === "Delete" || e.key === "Backspace") {
    e.preventDefault();
    if (state.selection.length) serverEdit({ type: "delete", atoms: state.selection });
    else deleteTarget(state.target);
    return;
  }
  if (e.key === "Escape") {
    if (drag) {
      drag = null;
      viewer.setGhost(null);
      $("selbox").hidden = true;
    } else setSelection([]);
    return;
  }
  if (hotkeyApplies(e.key)) {
    e.preventDefault();
    if (!e.repeat) hotkey(e.key);
    return;
  }
  if (e.key === "1") setTool("select");
  else if (e.key === "2") setTool("draw");
  else if (e.key === "3") setTool("delete");
  else if (key === "r") viewer.fit(model(), false);
  else if (e.key === "+" || e.key === "=") changeCharge(state.selection, +1);
  else if (e.key === "-") changeCharge(state.selection, -1);
  else if (KEY_ELEMENTS[e.key]) setElement(KEY_ELEMENTS[e.key]);
});

window.addEventListener("beforeunload", (e) => {
  if (hasUnsavedChanges()) {
    e.preventDefault();
    e.returnValue = "";
  }
});

// ---------------------------------------------------------- toolbar

for (const el of QUICK_ELEMENTS) {
  const b = document.createElement("button");
  b.textContent = el;
  b.dataset.element = el;
  b.title = `Draw ${el}`;
  b.addEventListener("click", () => setElement(el));
  $("elements").appendChild(b);
}
$("element-other").addEventListener("change", (e) => {
  const value = e.target.value.trim();
  if (!value) return;
  const sym = normalizeSymbol(value);
  if (!sym) {
    setMessage(`"${value}" is not an element symbol.`, "error");
    return;
  }
  e.target.value = sym;
  e.target.blur();
  setElement(sym);
  setMessage(`Drawing ${sym}.`);
});
for (const b of document.querySelectorAll("[data-tool]")) b.addEventListener("click", () => setTool(b.dataset.tool));

$("mode-3d").addEventListener("click", () => switchMode(3));
$("mode-2d").addEventListener("click", () => switchMode(2));
$("btn-open").addEventListener("click", () => canEdit() && $("file-input").click());
$("file-input").addEventListener("change", (e) => {
  const file = e.target.files[0];
  e.target.value = ""; // allow re-opening the same file
  openFile(file);
});
$("btn-save").addEventListener("click", save);
$("btn-undo").addEventListener("click", undo);
$("btn-redo").addEventListener("click", redo);
$("smiles-form").addEventListener("submit", (e) => {
  e.preventDefault();
  insertSmiles($("smiles").value);
});
$("btn-addh").addEventListener("click", addHydrogens);
$("btn-opt").addEventListener("click", optimize);
for (const id of ["method", "save-format"]) {
  $(id).addEventListener("change", (e) => {
    updateUI();
    e.target.blur(); // hand the keyboard back to the view (hotkeys)
  });
}
$("btn-reset").addEventListener("click", () => viewer.fit(model(), false));
$("btn-clear").addEventListener("click", () => {
  if (!canEdit()) return;
  state.selection = [];
  commit(withDim(M.emptyModel(), dim()), { view: "reset" });
});
$("btn-charge-minus").addEventListener("click", () => changeCharge(state.selection, -1));
$("btn-charge-plus").addEventListener("click", () => changeCharge(state.selection, +1));
$("keep-h").addEventListener("change", (e) => {
  state.keepH = e.target.checked;
  saveSetting("keepH", state.keepH);
  e.target.blur();
});
$("hotkeys").addEventListener("change", (e) => {
  state.hotkeys = e.target.checked;
  saveSetting("hotkeys", state.hotkeys);
  e.target.blur();
});
// Buttons keep focus after a click; blur them so Space/Enter do not re-trigger them.
document.addEventListener("click", (e) => {
  if (e.target instanceof HTMLButtonElement) e.target.blur();
});

function keyList(box, items) {
  box.textContent = "";
  // Keys that share a name are listed together, e.g. "3 a phenyl".
  const groups = new Map();
  for (const it of items) {
    if (!groups.has(it.name)) groups.set(it.name, []);
    groups.get(it.name).push(it.key);
  }
  for (const [name, keys] of groups) {
    const span = document.createElement("span");
    for (const k of keys) {
      const kbd = document.createElement("kbd");
      kbd.textContent = k;
      span.appendChild(kbd);
    }
    span.appendChild(document.createTextNode(` ${name}`));
    box.appendChild(span);
  }
}

function showKeys(fragments, bondKeys) {
  for (const f of fragments) state.fragments.set(f.key, f);
  for (const b of bondKeys) state.bondKeys.set(b.key, b);
  keyList($("keys-attach"), fragments.filter((f) => f.mode === "attach"));
  keyList($("keys-replace"), fragments.filter((f) => f.mode === "replace"));
  keyList($("keys-bond"), bondKeys);
}

// ------------------------------------------------------------ start-up

function setMethods(info, select = null) {
  const menu = $("method");
  const keep = select || menu.value;
  menu.textContent = "";
  for (const m of info.methods) menu.add(new Option(m, m));
  if (!info.methods.includes("g-xTB")) {
    const opt = new Option("g-xTB (not installed)", "");
    opt.disabled = true;
    opt.title = info.gxtb_installable
      ? "Use “Install g-xTB…” next to this menu."
      : "Install g-xTB from github.com/grimme-lab/g-xtb as 'gxtb' on PATH, or set MINIMOL_GXTB.";
    menu.add(opt);
  }
  if (info.methods.includes(keep)) menu.value = keep;
  $("btn-gxtb").hidden = !info.gxtb_installable;
  updateUI();
}

async function installGxtb() {
  if (!canEdit()) return;
  const ok = window.confirm(
    "minimol will download g-xTB 2.0.1 (about 40 MB) from github.com/grimme-lab/g-xtb, " +
      "check it against its published checksum and install it for your user account.\n\n" +
      "g-xTB is GPL-3.0 software by the Grimme group (University of Bonn); it is not part of minimol. Continue?",
  );
  if (!ok) return;
  const r = await busyCall("Downloading and installing g-xTB…", () => api("/api/install_gxtb"));
  if (!r) return;
  setMethods(r, "g-xTB");
  setMessage(`g-xTB installed (${r.path}).`);
}
$("btn-gxtb").addEventListener("click", installGxtb);

async function start() {
  updateUI();
  let info;
  try {
    info = await api("/api/info");
  } catch (err) {
    setMessage(err.message, "error");
    return;
  }
  showKeys(info.fragments || [], info.bond_keys || []);
  setMethods(info);
  try {
    const r = await api("/api/initial");
    if (r.molecule && commit(r.molecule, { view: "reset" })) {
      state.fileName = r.name.replace(/\.[^.]*$/, "") || "molecule";
      state.saved = model();
      updateUI();
      showNotes(r.notes, `Opened ${r.name}.`);
    }
  } catch (err) {
    setMessage(err.message, "error");
  }
}

// Used by the desktop window to ask before closing.
window.minimolHasUnsavedChanges = hasUnsavedChanges;
// Exposed for automated tests only.
window.__minimol = {
  state,
  history,
  viewers,
  get viewer() {
    return viewer;
  },
  model,
  commit,
  M,
};

start();
