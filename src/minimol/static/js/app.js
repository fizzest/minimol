// minimol editor: wires the toolbar, mouse tools and keyboard to the model,
// the viewer and the server.
//
// Rules that keep the editor predictable:
//  * history.current is the only molecule; every change goes through commit().
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

const $ = (id) => document.getElementById(id);
const CLICK_SLOP = 4; // pixels a pointer may move and still count as a click
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

let viewer;
try {
  viewer = new Viewer($("viewport"));
} catch (err) {
  window.minimolWebGLFailed = true;
  $("hint").textContent =
    `3D graphics (WebGL) could not start: ${err.message} ` +
    "In the desktop app, try starting it with: minimol --software-gl";
  throw err;
}
viewer.setModel(history.current);
viewer.fit(history.current);

const model = () => history.current;

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

function refreshModel(view = null) {
  const m = model();
  state.selection = state.selection.filter((i) => i < m.atoms.length);
  viewer.setModel(m);
  viewer.setSelection(state.selection);
  if (view === "reset") viewer.fit(m, false);
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
  const body = { molecule: m, op, keep_h: state.keepH, view: viewer.viewAxes() };
  const r = await busyCall("Editing…", () => api("/api/edit", body), { quiet: true });
  if (!r || model() !== m) return null;
  state.selection = state.selection.map((i) => r.map[i]).filter((i) => i >= 0);
  if (!commit(r.molecule)) return null;
  if (r.notes.length) setMessage(r.notes.join(" "), "warn");
  else setMessage("");
  return r;
}

async function openFile(file) {
  if (!file || !canEdit()) return;
  if (file.size > 50 * 1024 * 1024) {
    setMessage("That file is too large (over 50 MB).", "error");
    return;
  }
  const text = await file.text();
  const r = await busyCall(`Opening ${file.name}…`, () => api("/api/import", { name: file.name, text }));
  if (!r) return;
  state.selection = [];
  if (commit(r.molecule, { view: "reset" })) {
    state.fileName = file.name.replace(/\.[^.]*$/, "") || "molecule";
    state.saved = model();
    updateUI();
    showNotes(r.notes, `Opened ${file.name}.`);
  }
}

async function save() {
  const m = model();
  if (m.atoms.length === 0) {
    setMessage("Nothing to save yet.", "warn");
    return;
  }
  const fmt = $("save-format").value;
  const r = await busyCall("Saving…", () => api("/api/export", { molecule: m, format: fmt, title: state.fileName }));
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

let pendingSave = null;

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
  const r = await busyCall("Building 3D structure…", () => api("/api/smiles", { smiles }));
  if (!r) return;
  const m = model();
  if (m.atoms.length === 0) {
    commit(r.molecule, { view: "reset" });
  } else {
    commit(M.append(m, r.molecule, M.besideOffset(m, r.molecule)).model, { view: "fit" });
  }
  $("smiles").value = "";
  showNotes(r.notes, `Inserted ${smiles}.`);
}

async function addHydrogens() {
  if (!canEdit()) return;
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
  if (!t) return false;
  if (t.kind === "bond") {
    if (!["1", "2", "3"].includes(key)) return false;
    const b = model().bonds[t.index];
    await serverEdit({ type: "bond_order", a: b.a, b: b.b, order: Number(key) });
    return true;
  }
  if (key === "+" || key === "=") {
    changeCharge([t.index], +1);
    return true;
  }
  if (key === "-") {
    changeCharge([t.index], -1);
    return true;
  }
  if (!state.fragments.has(key)) return false;
  const r = await serverEdit({ type: "fragment", atom: t.index, key });
  if (r && r.focus !== null) setTarget({ kind: "atom", index: r.focus });
  return true;
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
    return `${a.el}${pick.index + 1} · ${nb} bond${nb === 1 ? "" : "s"}${a.charge ? ` · charge ${a.charge > 0 ? "+" : ""}${a.charge}` : ""}`;
  }
  const b = m.bonds[pick.index];
  return `bond ${atomLabel(b.a)}–${atomLabel(b.b)} · order ${b.order} · ${M.distance(M.pos(m.atoms[b.a]), M.pos(m.atoms[b.b])).toFixed(3)} Å`;
}

function measurementText() {
  const m = model();
  const s = state.selection;
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

  const locked = state.busy && !$("busy").hidden;
  for (const el of document.querySelectorAll("[data-edit], #method, #uhf, #fix-selected, #smiles, #element-other")) {
    el.disabled = locked;
  }
  if (!locked) {
    $("btn-undo").disabled = !history.canUndo;
    $("btn-redo").disabled = !history.canRedo;
  }
  $("uhf-wrap").hidden = $("method").value !== "g-xTB";
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
const canvas = viewer.canvas;

canvas.addEventListener("pointerdown", (e) => {
  if (drag) return;
  canvas.setPointerCapture(e.pointerId);
  const pick = viewer.pick(e.clientX, e.clientY);
  let mode = "rotate";
  if (e.button === 0 && (e.ctrlKey || e.metaKey)) mode = "box";
  else if (e.button === 1 || e.button === 2 || e.shiftKey) mode = "pan";
  else if (state.tool === "draw" && pick && pick.kind === "atom" && !state.busy) mode = "bond";
  drag = { id: e.pointerId, button: e.button, x0: e.clientX, y0: e.clientY, x: e.clientX, y: e.clientY, mode, pick, moved: false };
});

canvas.addEventListener("pointermove", (e) => {
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
});

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
  if (canvas.hasPointerCapture(e.pointerId)) canvas.releasePointerCapture(e.pointerId);
  viewer.setGhost(null);
  $("selbox").hidden = true;
  if (cancelled) return;
  if (d.mode === "box") {
    if (d.moved) {
      const inside = viewer.atomsInRect(d.x0, d.y0, e.clientX, e.clientY);
      setSelection(e.shiftKey ? [...state.selection, ...inside] : inside);
    } else if (d.pick) {
      const m = model();
      toggleSelected(d.pick.kind === "atom" ? [d.pick.index] : [m.bonds[d.pick.index].a, m.bonds[d.pick.index].b]);
    }
    return;
  }
  if (!d.moved && d.button === 0) handleClick(d.pick, e);
  else if (!d.moved && d.button === 2 && state.tool !== "select") deleteTarget(d.pick);
  else if (d.moved && d.mode === "bond") handleBondDrag(d.pick.index, e);
}

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
    viewer.zoom(Math.exp(Math.max(-200, Math.min(200, e.deltaY * scale)) * 0.0015));
  },
  { passive: false },
);

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

window.addEventListener("keydown", (e) => {
  const t = e.target;
  const typing = t instanceof HTMLElement && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT") && t.type !== "checkbox";
  const mod = e.ctrlKey || e.metaKey;
  const key = e.key.toLowerCase();
  if (typing) {
    if (e.key === "Escape") t.blur();
    return;
  }
  if (mod) {
    const actions = {
      z: () => (e.shiftKey ? redo() : undo()),
      y: redo,
      s: save,
      r: optimize,
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
  if (state.hotkeys && state.target && !e.repeat) {
    const known = state.target.kind === "bond" ? ["1", "2", "3"].includes(e.key) : state.fragments.has(e.key) || "+=-".includes(e.key);
    if (known) {
      e.preventDefault();
      hotkey(e.key);
      return;
    }
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
  setElement(sym);
  setMessage(`Drawing ${sym}.`);
});
for (const b of document.querySelectorAll("[data-tool]")) b.addEventListener("click", () => setTool(b.dataset.tool));

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
$("method").addEventListener("change", updateUI);
$("btn-reset").addEventListener("click", () => viewer.fit(model(), false));
$("btn-clear").addEventListener("click", () => {
  if (!canEdit()) return;
  state.selection = [];
  commit(M.emptyModel(), { view: "reset" });
});
$("btn-charge-minus").addEventListener("click", () => changeCharge(state.selection, -1));
$("btn-charge-plus").addEventListener("click", () => changeCharge(state.selection, +1));
$("keep-h").addEventListener("change", (e) => {
  state.keepH = e.target.checked;
  saveSetting("keepH", state.keepH);
  e.target.blur(); // give keys back to the view
});
$("hotkeys").addEventListener("change", (e) => {
  state.hotkeys = e.target.checked;
  saveSetting("hotkeys", state.hotkeys);
  e.target.blur();
});

function showFragmentKeys(fragments) {
  for (const f of fragments) state.fragments.set(f.key, f);
  for (const mode of ["attach", "replace"]) {
    const box = $(`keys-${mode}`);
    box.textContent = "";
    // Keys that share a name are listed together, e.g. "3 a phenyl".
    const groups = new Map();
    for (const f of fragments.filter((x) => x.mode === mode)) {
      if (!groups.has(f.name)) groups.set(f.name, []);
      groups.get(f.name).push(f.key);
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
}

// ------------------------------------------------------------ start-up

async function start() {
  updateUI();
  let info;
  try {
    info = await api("/api/info");
  } catch (err) {
    setMessage(err.message, "error");
    return;
  }
  showFragmentKeys(info.fragments || []);
  const select = $("method");
  for (const m of info.methods) select.add(new Option(m, m));
  if (!info.methods.includes("g-xTB")) {
    const opt = new Option("g-xTB (not installed)", "");
    opt.disabled = true;
    opt.title = "Install g-xTB from github.com/grimme-lab/g-xtb as 'gxtb' on PATH, or set MINIMOL_GXTB.";
    select.add(opt);
  }
  updateUI();
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
window.__minimol = { state, history, viewer, model, commit, M };

start();
