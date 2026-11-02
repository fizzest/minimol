// minimol editor: wires the toolbar, mouse tools and keyboard to the model,
// the viewer and the server.
//
// Rules that keep the editor predictable:
//  * history.current is the only molecule; every change goes through commit().
//  * While a server call runs (state.busy) nothing may change the molecule, so
//    a result always applies to the structure it was computed for.
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
  hover: null,
  busy: false,
  energy: null, // { model, text }
  fileName: "molecule",
  saved: null, // model last saved or opened
};

const history = new History(M.emptyModel());
state.saved = history.current;
const viewer = new Viewer($("viewport"));
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
  state.hover = null;
  $("hover").textContent = "";
  viewer.setModel(m);
  viewer.setHover(null);
  viewer.setSelection(state.selection);
  if (view === "reset") viewer.fit(m, false);
  else if (view === "fit") viewer.fit(m, true);
  updateUI();
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

function deleteAtoms(indices) {
  if (!canEdit() || indices.length === 0) return;
  const { model: next, map } = M.removeAtomsWithHydrogens(model(), indices);
  state.selection = state.selection.map((i) => map[i]).filter((i) => i >= 0);
  commit(next);
}

function changeCharge(delta) {
  if (!canEdit()) return;
  if (state.selection.length === 0) {
    setMessage("Select atoms first (Select tool), then change their charge.", "warn");
    return;
  }
  commit(M.changeCharge(model(), state.selection, delta));
}

// ------------------------------------------------------------- display

function atomLabel(i) {
  return `${model().atoms[i].el}${i + 1}`;
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

function updateUI() {
  const m = model();
  const n = m.atoms.length;
  const q = M.totalCharge(m);
  $("info").textContent = n ? `${M.formula(m)} · ${n} atom${n === 1 ? "" : "s"} · charge ${q > 0 ? "+" : ""}${q}` : "Empty";
  $("measure").textContent = measurementText();
  $("energy").textContent = state.energy && state.energy.model === m ? state.energy.text : "";
  $("hint").hidden = n > 0;

  for (const b of document.querySelectorAll("[data-tool]")) b.classList.toggle("active", b.dataset.tool === state.tool);
  for (const b of document.querySelectorAll("[data-element]")) b.classList.toggle("active", b.dataset.element === state.element);
  $("element-other").classList.toggle("active", !QUICK_ELEMENTS.includes(state.element));
  $("viewport").className = `tool-${state.tool}`;

  for (const el of document.querySelectorAll("[data-edit], #method, #uhf, #fix-selected, #smiles, #element-other")) {
    el.disabled = state.busy;
  }
  if (!state.busy) {
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
  if (state.tool !== "draw") state.tool = "draw";
  if (QUICK_ELEMENTS.includes(el)) $("element-other").value = "";
  updateUI();
}

// --------------------------------------------------- server round trips

async function busyCall(text, fn) {
  if (state.busy) return null;
  state.busy = true;
  $("busy-text").textContent = text;
  $("busy").hidden = false;
  updateUI();
  try {
    return await fn();
  } catch (err) {
    setMessage(err.message, "error");
    return null;
  } finally {
    state.busy = false;
    $("busy").hidden = true;
    updateUI();
  }
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
  state.saved = m;
  setMessage(`Saved ${a.download}.`);
}

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
  if (!r) return;
  if (model() !== m) return;
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

// ------------------------------------------------------- mouse tools

function hoverText(pick) {
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
    const additive = e.shiftKey || e.ctrlKey || e.metaKey;
    if (!pick) {
      if (!additive) setSelection([]);
    } else if (pick.kind === "atom") {
      const i = pick.index;
      if (additive) {
        setSelection(state.selection.includes(i) ? state.selection.filter((j) => j !== i) : [...state.selection, i]);
      } else {
        setSelection([i]);
      }
    } else {
      const b = m.bonds[pick.index];
      setSelection(additive ? [...state.selection, b.a, b.b] : [b.a, b.b]);
    }
    return;
  }
  if (!canEdit()) return;
  if (state.tool === "draw") {
    if (!pick) {
      commit(M.addAtom(m, state.element, viewer.screenToModel(e.clientX, e.clientY, viewer.centerPoint())).model);
    } else if (pick.kind === "atom") {
      const i = pick.index;
      if (m.atoms[i].el !== state.element) {
        commit(M.setElement(m, [i], state.element));
      } else {
        const p = M.newAtomPosition(m, i, state.element, viewer.viewAxes());
        const { model: next, index } = M.addAtom(m, state.element, p);
        commit(M.addBond(next, i, index, 1));
      }
    } else {
      commit(M.cycleBondOrder(m, pick.index));
    }
  } else if (state.tool === "delete") {
    if (pick && pick.kind === "atom") deleteAtoms([pick.index]);
    else if (pick) commit(M.removeBond(m, pick.index));
  }
}

function handleBondDrag(i, e) {
  if (!canEdit()) return;
  const m = model();
  const target = viewer.pick(e.clientX, e.clientY);
  if (target && target.kind === "atom") {
    const j = target.index;
    if (j === i) return;
    const k = M.findBond(m, i, j);
    commit(k >= 0 ? M.cycleBondOrder(m, k) : M.addBond(m, i, j, 1));
    return;
  }
  const p = dragTarget(i, e.clientX, e.clientY);
  if (!p) return;
  const { model: next, index } = M.addAtom(m, state.element, p);
  commit(M.addBond(next, i, index, 1));
}

let drag = null;
const canvas = viewer.canvas;

canvas.addEventListener("pointerdown", (e) => {
  if (drag) return;
  canvas.setPointerCapture(e.pointerId);
  const pick = viewer.pick(e.clientX, e.clientY);
  let mode = "rotate";
  if (e.button === 1 || e.button === 2 || e.shiftKey) mode = "pan";
  else if (state.tool === "draw" && pick && pick.kind === "atom" && !state.busy) mode = "bond";
  drag = { id: e.pointerId, button: e.button, x0: e.clientX, y0: e.clientY, x: e.clientX, y: e.clientY, mode, pick, moved: false };
});

canvas.addEventListener("pointermove", (e) => {
  if (!drag) {
    const pick = viewer.pick(e.clientX, e.clientY);
    state.hover = pick;
    viewer.setHover(pick);
    $("hover").textContent = hoverText(pick);
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

function endDrag(e, cancelled) {
  if (!drag || e.pointerId !== drag.id) return;
  const d = drag;
  drag = null;
  if (canvas.hasPointerCapture(e.pointerId)) canvas.releasePointerCapture(e.pointerId);
  viewer.setGhost(null);
  if (cancelled) return;
  if (!d.moved && d.button === 0) handleClick(d.pick, e);
  else if (d.moved && d.mode === "bond") handleBondDrag(d.pick.index, e);
}

canvas.addEventListener("pointerup", (e) => endDrag(e, false));
canvas.addEventListener("pointercancel", (e) => endDrag(e, true));
canvas.addEventListener("pointerleave", () => {
  if (!drag) {
    state.hover = null;
    viewer.setHover(null);
    $("hover").textContent = "";
  }
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
  if (mod && key === "z") {
    e.preventDefault();
    if (e.shiftKey) redo();
    else undo();
  } else if (mod && key === "y") {
    e.preventDefault();
    redo();
  } else if (mod && key === "s") {
    e.preventDefault();
    save();
  } else if (mod && key === "o") {
    e.preventDefault();
    if (canEdit()) $("file-input").click();
  } else if (mod && key === "a") {
    e.preventDefault();
    setSelection(model().atoms.map((_, i) => i));
  } else if (mod || e.altKey) {
    return;
  } else if (e.key === "Delete" || e.key === "Backspace") {
    e.preventDefault();
    deleteAtoms(state.selection);
  } else if (e.key === "Escape") {
    if (drag) {
      drag = null;
      viewer.setGhost(null);
    } else setSelection([]);
  } else if (e.key === "1") setTool("select");
  else if (e.key === "2") setTool("draw");
  else if (e.key === "3") setTool("delete");
  else if (key === "r") viewer.fit(model(), false);
  else if (e.key === "+" || e.key === "=") changeCharge(+1);
  else if (e.key === "-") changeCharge(-1);
  else if (KEY_ELEMENTS[key]) setElement(KEY_ELEMENTS[key]);
});

window.addEventListener("beforeunload", (e) => {
  if (model() !== state.saved && model().atoms.length > 0) {
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
$("btn-charge-minus").addEventListener("click", () => changeCharge(-1));
$("btn-charge-plus").addEventListener("click", () => changeCharge(+1));

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
      showNotes(r.notes, `Opened ${r.name}.`);
    }
  } catch (err) {
    setMessage(err.message, "error");
  }
}

// Exposed for automated tests only.
window.__minimol = { state, history, viewer, model, commit, M };

start();
