# minimol

A deliberately small 3D molecule viewer and editor for preparing structures:
look at molecules, build or fix them in 3D (with ChemDraw-style fragment
hotkeys), optimize their geometry with **MMFF94**, **MMFF94s**, **UFF** or,
optionally, **g-xTB**, and save them for your quantum-chemistry program.

It runs as a desktop application (its own window). Nothing is sent anywhere:
the window talks only to the minimol process on your own machine.

## Install and run

Requires Python 3.10 or newer.

```bash
pip install .
```

This also installs RDKit and PySide6 (Qt). PySide6 is large: about 670 MB
once installed.

```bash
minimol                    # open the window
```

```bash
minimol molecule.xyz       # open a file straight away
```

* `--browser` uses your web browser instead of the window. Ctrl+R optimizes
  there too in Chromium (tested); other browsers may reload the page instead,
  so use the Optimize button.
* `--software-gl` lets the window draw without a usable GPU (virtual
  machines, remote desktops). Use it if the window says 3D graphics could not
  start. This was verified on Linux only.

### Optional: g-xTB

[g-xTB](https://github.com/grimme-lab/g-xtb) (Grimme group, preprint:
Froitzheim, Müller, Hansen, Grimme, ChemRxiv 2025, DOI
10.26434/chemrxiv-2025-bjxv) is a semiempirical method for H–Lr. It is
distributed as a modified `xtb` binary under the GPL-3.0. minimol does **not**
bundle it; it runs it as a separate program if you install it:

1. Download the archive for your OS from the `binaries/` folder of
   https://github.com/grimme-lab/g-xtb and check it against its `.sha256` file.
2. Extract it and make its `bin/xtb` reachable as `gxtb`, for example on Linux:

   ```bash
   mkdir -p ~/.local/opt ~/.local/bin
   tar xJf xtb-6.7.1-gxtb-*-linux-x86_64.tar.xz -C ~/.local/opt
   ln -s ~/.local/opt/xtb-6.7.1/bin/xtb ~/.local/bin/gxtb
   ```

   Alternatively set `MINIMOL_GXTB=/full/path/to/xtb` (the way to do it on
   Windows). On macOS and Windows the bundled shared libraries must stay next
   to the binary (see the g-xTB README).
3. Restart minimol. "g-xTB" now appears in the method list.

minimol only accepts a binary whose `--help` lists `--gxtb`. It runs
`xtb in.xyz --gxtb --opt --norestart --chrg <sum of formal charges>` in a
temporary folder; atoms marked "fix selected" become an xcontrol `$fix` block.
Tested with g-xTB 2.0.1 (xtb 6.7.1, built 2026-05-14) on Linux x86-64.

## Using it

### Mouse

| Action | How |
| --- | --- |
| Rotate / pan / zoom | left-drag on empty space (3D) / right-drag or Shift+drag / wheel or pinch |
| New atom | **Draw** tool: click empty space |
| Grow a bond | Draw: click an atom, or drag from it (the new atom goes where the preview line shows) |
| Bond two atoms | Draw: drag from one atom onto another |
| Change element | pick an element (or type one, e.g. `Fe`), then click an atom |
| Bond order | Draw: click a bond to cycle single → double → triple |
| Delete | **right-click** an atom or bond (Draw or Delete tool), or click it with the Delete tool |
| Select | **Ctrl+drag** a box (in 2D also a plain drag on empty space; Ctrl+Shift+drag adds), Ctrl+click toggles an atom; Select tool: click, Shift+click adds |
| Measure | select 2, 3 or 4 atoms for distance, angle, dihedral |

### Keys

| Key | Action |
| --- | --- |
| **Ctrl+R** | optimize the geometry (3D) / clean up the drawing (2D) |
| **Tab** | switch between 3D and 2D |
| Ctrl+C / Ctrl+X / Ctrl+V | copy, cut, paste structures (see Clipboard) |
| Ctrl+Z / Ctrl+Shift+Z (Ctrl+Y) | undo / redo |
| Delete | delete the selection, or the atom/bond under the mouse |
| Ctrl+A, Esc | select all, clear selection |
| Ctrl+O, Ctrl+S | open, save |
| 1 / 2 / 3, R | Select / Draw / Delete tool, reset view (when not hovering an atom) |

On macOS use Cmd instead of Ctrl.

### Fragment hotkeys (ChemDraw-style)

Hover an atom or bond and press a key; this works in 3D and 2D. The key map
follows ChemDraw's hotkey cheat sheet, but it is minimol's own and may differ
in details. The sidebar lists every key.

* **On an atom, most keys turn that atom into a group**, as in ChemDraw: `3`
  on a CH₃ makes it a phenyl ring, `o` makes it OH.
  * Rings and branches: `3`/`a`/`P` phenyl, `6` cyclohexyl, `7` cyclopentyl,
    `v` cyclopropyl, `u` cyclobutyl, `9` isopropyl, `K` tert-butyl.
  * Unsaturated groups: `2` C=O, `8` vinyl, `z` ethynyl, `Z` azide, `k` SO₂.
  * Labels: `o`/`q` OH, `O` OMe, `n`/`w` NH₂, `N` NO₂, `s` SH, `S` SiH₃,
    `p` PH₂, `f` F, `F` CF₃, `l`/`C` Cl, `L` Li, `b` Br, `B` BH₂, `i` I, `h` H,
    `m`/`c` CH₃, `e` Et, `A` Ac, `E` CO₂Me, `M` MgBr, `y` Boc, `H` Cbz,
    `Q` Fmoc.
  * Where the atom goes:
    * On a hydrogen, the hydrogen becomes the group.
    * On an atom inside a chain or ring, single-atom groups change the
      element.
    * Larger groups are added there as a substituent instead (you get a
      note).
* **Keys that add a new atom to the hovered one:** `1` carbon, `4` methyl
  towards you, `5` methyl away from you (the 3D meaning of ChemDraw's
  wedge/hash keys; in 2D both just add a methyl).
* **On a bond:**
  * `1`/`2`/`3` set the order.
  * `a` fuses a benzene ring.
  * `v`, `4`, `5`, `6`, `7`, `8` fuse a 3- to 8-membered ring.
  * Rings go on the side of the bond with the most room.
* **On an atom:** `+`/`−` change the charge, Delete removes it.

**The target moves on after each key**, so keys chain:
* On methane, `1 1 o` adds two carbons and turns the last into OH, giving
  ethanol.
* After `3` the target is the ring's para carbon.
* Moving the mouse picks a new target.

Groups are placed at standard bond lengths and turned to avoid clashes, but
not optimized: press Ctrl+R when done. The "Fragment hotkeys" checkbox turns
hotkeys off.

### 2D mode (ChemDraw-style drawing)

The **3D / 2D** switch in the toolbar (or Tab) toggles between the 3D editor
and a skeletal drawing:

* **Drawing style:** carbons are plain vertices and hydrogens are implicit,
  shown in labels such as OH or NH₂. Atoms with too many bonds get a red box.
* **The same tools and hotkeys** work in 2D. Bonds have one length and snap to
  30° steps when dragged, and chains zigzag at 120°.
* **Ctrl+R** cleans the drawing up (RDKit's 2D layout).
* **Going to 3D** adds hydrogens, embeds the structure (ETKDG) and relaxes it
  with MMFF94, or UFF if MMFF94 lacks parameters. A drawing has no wedges, so
  stereocentres come out arbitrary and you get a note. **Going back to 2D**
  lays the structure out again.
* **Switching straight back without editing** restores the exact previous
  drawing or 3D structure.
* **Undo/redo** works across switches.

### Hydrogens

* **Bonding to a full atom replaces one of its hydrogens:** the one pointing
  where the new bond goes. An atom with a free valence keeps its hydrogens.
* **Keep hydrogens (on by default):** new atoms come with hydrogens, and
  changing an element, bond order or charge, or deleting a neighbour,
  re-derives the hydrogens of the atoms involved (RDKit's valence model).
  Hydrogens elsewhere are never touched. Turn it off to draw bare skeletons;
  **Add H** then fills every open valence at once.

**Energies** are what each method reports: kcal/mol for the force fields,
Hartree (Eh) for g-xTB. The energy display disappears as soon as the
structure changes, so it always belongs to what is on screen. The total
charge is the sum of formal charges. For g-xTB you can set the number of
unpaired electrons; blank means singlet for an even electron count and
doublet for an odd one.

### Clipboard (copy/paste with ChemDraw)

* **Copy** (Ctrl+C) puts the selection, or the whole molecule, on the
  clipboard in several formats:
  * SMILES as plain text
  * a 2D molfile (`chemical/x-mdl-molfile`)
  * CDXML (`chemical/x-cdxml`)
  * MDL CT under the Windows clipboard format name `MDLCT`

  Hydrogens are implicit and the coordinates are a 2D drawing. In ChemDraw,
  *Paste Special › SMILES* reads the text. Whether ChemDraw's plain Ctrl+V
  picks up one of the other formats depends on your ChemDraw version and
  settings. **This has not been tested against ChemDraw.**
* **Paste** (Ctrl+V) reads whatever readable format is on the clipboard,
  checking the data itself rather than the format names:
  * CDXML
  * MDL CT (what ChemDraw and ISIS/Draw put on the Windows clipboard)
  * molfile text
  * SMILES
  * InChI

  The structure is added next to what is on screen: as a drawing in 2D, or
  built in 3D in 3D mode.
* **Binary CDX** (ChemDraw's native format) is recognised. **The RDKit
  release this was built with (2026.09) cannot read it from Python,** so if
  that is all the clipboard holds you get a message listing the formats that
  were there:
  * In ChemDraw, copy as SMILES or MOL text instead (Edit › Copy As), or turn
    on MDL CT clipboard output in its preferences.
  * minimol uses RDKit's CDX reader automatically once an RDKit version
    provides it.
* **In browser mode** only plain text (SMILES, molfile, InChI) goes through
  the clipboard.

### Files

* **Open:** `.xyz`, `.mol`, `.sdf`, `.mol2`, `.pdb`, `.cdxml`, via the button,
  Ctrl+O or drag and drop.
  * Only the first molecule or frame of a multi-record file is loaded, and you
    are told when that happens; all structures in a CDXML file are loaded.
  * In 2D mode, files with 2D coordinates keep their layout.
  * Binary `.cdx` files cannot be opened (see Clipboard); save as `.cdxml` in
    ChemDraw.
* **XYZ bonds** are perceived with RDKit's `DetermineBonds`, assuming a
  neutral molecule. If that fails, all bonds become single and you get a note.
  Above 200 atoms only connectivity is perceived.
* **2D MOL/SDF files** are turned into 3D (hydrogens added, relaxed with
  MMFF94 or UFF).
* **Save:** `.xyz` (8 decimals), `.sdf`/`.mol` (4 decimals, a limit of the
  format), `.pdb` or `.cdxml`, through a native save dialog.
  * 2D drawings save as `.mol`, `.sdf` or `.cdxml` with implicit hydrogens.
  * A 3D structure saved as `.cdxml` is laid out in 2D.
* **SMILES:** type one in the toolbar and press Enter.
  * In 3D it is embedded (RDKit ETKDG, fixed seed) and relaxed.
  * In 2D it is laid out flat.
  * Either way it is placed next to anything already there.

## How it is built

* **Chemistry:** [RDKit](https://www.rdkit.org) (BSD-3-Clause): force fields,
  3D embedding, valences and hydrogens, bond perception, file formats.
  * `src/minimol/edit.py`: every structural edit, fragment and ring fusion,
    in 3D and 2D.
  * `src/minimol/layout.py`: 2D drawings and 2D ↔ 3D conversion.
  * `src/minimol/clipboard.py`: clipboard formats (MDL CT, CDXML, ...).
  * `src/minimol/chem.py`: import, export and optimization.
  * Force-field optimizations always include interactions between separate
    fragments; RDKit's helpers default to leaving them out.
* **g-xTB:** an optional external program (`src/minimol/gxtb.py`).
* **Window:** PySide6/Qt WebEngine (`src/minimol/desktop.py`). It shows the
  editor page served by a local, stateless server (`server.py`, Python
  standard library, bound to 127.0.0.1, with Host and header checks). The
  window also bridges the server to the system clipboard.
* **Editor:** plain ES modules with no build step. The 3D view uses
  [three.js](https://threejs.org) r187 (MIT, vendored unmodified); the 2D view
  is a canvas drawing (`static/js/viewer2d.js`). The molecule is an immutable value
  and undo/redo is a list of them. Edits go to the server one at a time and
  editing is locked while one runs, so a result never lands on a structure
  that changed meanwhile.

## Development

```bash
python -m venv .venv && .venv/bin/pip install -e ".[test]"
```

```bash
.venv/bin/python -m pytest      # chemistry, edits, g-xTB, server, desktop window
```

```bash
npm test                        # model/undo logic (Node 22+)
```

```bash
npm install --no-save playwright && npx playwright install chromium
```

```bash
npm run e2e                     # drives the editor in Chromium with real mouse/keyboard input
```

* g-xTB tests skip when it is not installed. The desktop tests skip without
  PySide6 and run offscreen (`QT_QPA_PLATFORM=offscreen`).
* `scripts/gen_elements.py` regenerates `static/js/elements.js` from RDKit's
  periodic table.

## Not included (yet)

* Orbitals and densities (IboView-style isosurfaces).
* Stereo bonds (wedges and hashes) in 2D.
* Reading binary CDX (waiting for RDKit).
* Periodic systems.
* Conformer search.
* Frequencies.
* Solvent models.

## License

minimol is released under the MIT License (see [LICENSE](LICENSE)). It
includes three.js (MIT, `src/minimol/static/vendor/three-LICENSE.txt`) and
depends on RDKit (BSD-3-Clause) and PySide6 (LGPL-3.0). g-xTB (GPL-3.0) is
not included; minimol only runs it as a separate program if you install it.
