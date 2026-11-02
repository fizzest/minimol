# minimol

A deliberately small 3D molecule viewer and editor: look at molecules, build
or fix them in 3D, and optimize their geometry with **MMFF94**, **MMFF94s**,
**UFF**, or (optionally) **g-xTB**. Then save the structure for your
quantum-chemistry program.

It runs as a local program that opens in your web browser. Nothing is sent
anywhere: the page talks only to the minimol process on your own machine.

## Install and run

Requires Python 3.10 or newer.

```bash
pip install .
```

```bash
minimol                 # opens http://127.0.0.1:8765/ in your browser
```

```bash
minimol molecule.xyz    # open a file straight away
```

`minimol --port 9000` picks another port, and `--no-browser` skips opening a
browser. If the port is taken, a free one is used and printed.

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

   Alternatively set `MINIMOL_GXTB=/full/path/to/xtb`. On Windows, use the
   environment variable; on macOS and Windows the bundled shared libraries
   must stay next to the binary (see the g-xTB README).
3. Restart minimol. "g-xTB" now appears in the method list.

minimol only accepts a binary whose `--help` lists `--gxtb`, so a regular xtb
is never mistaken for g-xTB. The run uses `xtb in.xyz --gxtb --opt
--norestart --chrg <sum of formal charges>` in a temporary folder. Atoms
marked "fix selected" are passed as an xcontrol `$fix` block. g-xTB's
built-in optimizer does the minimization. This was tested with the g-xTB
2.0.1 binary (xtb 6.7.1, built 2026-05-14) on Linux x86-64.

## Using it

| Action | How |
| --- | --- |
| Rotate / pan / zoom | left-drag on empty space / right-drag or Shift+drag / wheel or pinch |
| New atom | **Draw** tool: click empty space |
| Grow a bond | Draw: click an atom (new atom points away from its neighbours), or drag from an atom |
| Bond two atoms | Draw: drag from one atom onto another |
| Change element | pick an element (or type one, e.g. `Fe`), then click an atom |
| Bond order | Draw: click a bond to cycle single → double → triple |
| Delete | **Delete** tool: click an atom (its hydrogens go too) or a bond; or select and press Del |
| Formal charge | **Select** atoms, then `+` / `−` (buttons or keys) |
| Measure | Select 2, 3 or 4 atoms (Shift+click) for distance, angle, dihedral |
| Hydrogens | **Add H** fills open valences; nothing is added behind your back |
| Optimize | choose a method, click **Optimize**; "fix selected" freezes the selected atoms |
| Undo / redo | Ctrl+Z / Ctrl+Shift+Z (or Ctrl+Y) |
| Other keys | 1/2/3 tools, H C N O F P S elements, R reset view, Ctrl+A select all, Ctrl+O open, Ctrl+S save, Esc clear selection |

**What you see is what you compute.** An atom with a missing bond is a radical
(you will get a note when optimizing). Click **Add H** first if you want a
closed-shell molecule. The total charge is the sum of the atoms' formal
charges; for g-xTB you can also set the number of unpaired electrons (blank
means singlet for an even electron count and doublet for an odd one).

**Energies** are what the method reports: kcal/mol for the force fields, in
Hartree (Eh) for g-xTB. The energy display disappears as soon as the
structure changes, so it always belongs to what is on screen.

### Files

* Open: `.xyz`, `.mol`, `.sdf`, `.mol2`, `.pdb` (button, Ctrl+O, or drag and
  drop onto the view). Only the first molecule or frame of a multi-record
  file is loaded, and you are told when that happens.
* XYZ files have no bonds, so they are perceived with RDKit's
  `DetermineBonds`, which assumes a neutral molecule. If bond orders cannot be
  assigned, all bonds become single and you get a note. Above 200 atoms only
  connectivity is perceived.
* MOL/SDF files with 2D coordinates are turned into 3D (hydrogens added,
  then relaxed with MMFF94, or UFF if MMFF94 lacks parameters).
* Save: `.xyz` (8 decimals), `.sdf`/`.mol` (4 decimals, a limit of the
  format), `.pdb`.
* SMILES: type one in the toolbar and press Enter. It is embedded in 3D
  (RDKit ETKDG, fixed random seed), relaxed, and placed next to anything
  already on screen.

## How it is built

* **Chemistry:** [RDKit](https://www.rdkit.org) (BSD-3-Clause): MMFF94/MMFF94s/
  UFF, 3D embedding, bond perception, file formats. Interactions between
  separate fragments are always included in force-field optimizations.
  RDKit's helpers default to leaving them out.
* **g-xTB:** optional external program, see above (`src/minimol/gxtb.py`).
* **Server:** Python standard library only (`src/minimol/server.py`). It
  listens on 127.0.0.1 only and is stateless: every request carries the
  whole molecule and gets a new one back. API calls need a custom header and
  a matching `Host`, so other web pages cannot use it.
* **Editor:** plain ES modules and [three.js](https://threejs.org) r187
  (MIT, vendored unmodified in `static/vendor/`), with no build step. The
  molecule is an immutable value: every edit in `static/js/model.js` is a
  pure function, and undo/redo is a list of those values. While a calculation
  runs, edits are refused, so a result can never land on a structure that
  changed in the meantime.

## Development

```bash
python -m venv .venv && .venv/bin/pip install -e ".[test]"
```

```bash
.venv/bin/python -m pytest          # chemistry, g-xTB, server (g-xTB tests skip if not installed)
```

```bash
npm test                            # model/undo logic (Node 22+)
```

```bash
npm install --no-save playwright && npx playwright install chromium
```

```bash
npm run e2e                         # drives the real editor in Chromium
```

`scripts/gen_elements.py` regenerates `static/js/elements.js` (symbols and
covalent radii) from RDKit's periodic table.

## Not included (yet)

Orbitals and densities (IboView-style isosurfaces), periodic systems,
conformer search, frequency calculations, and solvent models.
