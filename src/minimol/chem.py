"""Chemistry operations for minimol, implemented on top of RDKit.

Every function here is stateless: it takes a molecule in minimol's JSON form,
returns a new one, and never keeps anything between calls. The JSON form is

    {"atoms": [{"el": "C", "x": 0.0, "y": 0.0, "z": 0.0, "charge": 0}, ...],
     "bonds": [{"a": 0, "b": 1, "order": 1}, ...]}

Atom indices are positions in ``atoms``. Bond orders are 1, 2 or 3 (aromatic
rings are always written in a Kekulé form). Coordinates are in Ångström.

The editor is "what you see is what you get": a molecule built from JSON has
no implicit hydrogens, so an atom with an open valence is treated as a radical
by RDKit. Hydrogens are only added when the user asks for it
(:func:`add_hydrogens`).
"""

from __future__ import annotations

import math
import re

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds, rdDistGeom
from rdkit.Chem import rdForceFieldHelpers as FFH
from rdkit.Geometry import Point3D

# RDKit prints sanitization problems to stderr. We report them to the user
# through exceptions instead, so keep the console quiet.
RDLogger.DisableLog("rdApp.*")

FORCE_FIELDS = ("MMFF94", "MMFF94s", "UFF")
EXPORT_FORMATS = ("xyz", "mol", "sdf", "pdb")
IMPORT_EXTENSIONS = ("xyz", "mol", "sdf", "sd", "mol2", "pdb", "ent")

_ORDER_TO_TYPE = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
}
_TYPE_TO_ORDER = {v: k for k, v in _ORDER_TO_TYPE.items()}

# Bond-order perception from bare coordinates (xyz2mol) can take a very long
# time on large conjugated systems. Above this size only connectivity is
# perceived.
MAX_ATOMS_FOR_BOND_ORDERS = 200

_PT = Chem.GetPeriodicTable()


class ChemError(ValueError):
    """A problem with the user's molecule or input, safe to show to the user."""


# --------------------------------------------------------------------------
# JSON <-> RDKit
# --------------------------------------------------------------------------


def _atomic_number(symbol: object) -> int:
    if not isinstance(symbol, str) or not re.fullmatch(r"[A-Z][a-z]?[a-z]?", symbol):
        raise ChemError(f"Invalid element symbol: {symbol!r}")
    try:
        z = _PT.GetAtomicNumber(symbol)
    except RuntimeError:
        z = 0
    if z <= 0 or _PT.GetElementSymbol(z) != symbol:
        raise ChemError(f"Unknown element: {symbol!r}")
    return z


def _number(value: object, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ChemError(f"{what} must be a number")
    if not math.isfinite(value):
        raise ChemError(f"{what} must be finite")
    return float(value)


def _integer(value: object, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ChemError(f"{what} must be an integer")
    return value


def validate_json(data: object) -> None:
    """Raise :class:`ChemError` unless ``data`` is a well-formed molecule."""
    if not isinstance(data, dict):
        raise ChemError("Molecule must be an object")
    atoms = data.get("atoms")
    bonds = data.get("bonds", [])
    if not isinstance(atoms, list) or not isinstance(bonds, list):
        raise ChemError("Molecule needs 'atoms' and 'bonds' lists")
    for i, atom in enumerate(atoms):
        if not isinstance(atom, dict):
            raise ChemError(f"Atom {i} must be an object")
        _atomic_number(atom.get("el"))
        for key in ("x", "y", "z"):
            _number(atom.get(key), f"Atom {i} coordinate {key}")
        charge = _integer(atom.get("charge", 0), f"Atom {i} charge")
        if abs(charge) > 8:
            raise ChemError(f"Atom {i} charge {charge} is out of range")
    seen = set()
    n = len(atoms)
    for k, bond in enumerate(bonds):
        if not isinstance(bond, dict):
            raise ChemError(f"Bond {k} must be an object")
        a = _integer(bond.get("a"), f"Bond {k} atom a")
        b = _integer(bond.get("b"), f"Bond {k} atom b")
        order = _integer(bond.get("order", 1), f"Bond {k} order")
        if not (0 <= a < n and 0 <= b < n):
            raise ChemError(f"Bond {k} refers to a missing atom")
        if a == b:
            raise ChemError(f"Bond {k} connects atom {a} to itself")
        if order not in _ORDER_TO_TYPE:
            raise ChemError(f"Bond {k} has unsupported order {order}")
        key = (min(a, b), max(a, b))
        if key in seen:
            raise ChemError(f"Atoms {a} and {b} are bonded twice")
        seen.add(key)


def mol_from_json(data: object, *, sanitize: bool = True, implicit_h: bool = False) -> Chem.RWMol:
    """Build an RDKit molecule with one 3D conformer from minimol JSON.

    With ``implicit_h=False`` (the default) no atom gets implicit hydrogens,
    so the molecule is exactly what the user drew. ``implicit_h=True`` lets
    RDKit fill open valences of non-hydrogen atoms; it is used to work out
    how many hydrogens to add.
    """
    validate_json(data)
    atoms = data["atoms"]
    mol = Chem.RWMol()
    for atom in atoms:
        a = Chem.Atom(_atomic_number(atom["el"]))
        a.SetFormalCharge(atom.get("charge", 0))
        a.SetNoImplicit(not implicit_h or a.GetAtomicNum() == 1)
        mol.AddAtom(a)
    for bond in data.get("bonds", []):
        mol.AddBond(bond["a"], bond["b"], _ORDER_TO_TYPE[bond.get("order", 1)])
    conf = Chem.Conformer(len(atoms))
    for i, atom in enumerate(atoms):
        conf.SetAtomPosition(i, Point3D(float(atom["x"]), float(atom["y"]), float(atom["z"])))
    conf.Set3D(True)
    mol.AddConformer(conf, assignId=True)
    if sanitize:
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:  # RDKit raises several exception types here
            raise ChemError(_describe_sanitize_error(exc)) from None
    return mol


def _describe_sanitize_error(exc: Exception) -> str:
    # RDKit numbers atoms from 0; the UI (like xtb) numbers them from 1.
    text = str(exc).strip() or type(exc).__name__
    text = re.sub(r"atom # ?(\d+)", lambda m: f"atom {int(m.group(1)) + 1}", text)
    return f"Invalid structure: {text}"


def mol_to_json(mol: Chem.Mol, notes: list[str] | None = None) -> dict:
    """Convert an RDKit molecule (with a conformer) to minimol JSON."""
    if mol.GetNumConformers() == 0:
        raise ChemError("Molecule has no coordinates")
    work = Chem.Mol(mol)
    if any(b.GetIsAromatic() for b in work.GetBonds()):
        try:
            Chem.Kekulize(work, clearAromaticFlags=True)
        except Exception:
            raise ChemError("Could not assign a Kekulé structure to the aromatic rings") from None
    conf = work.GetConformer()
    atoms = []
    for atom in work.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        atoms.append(
            {
                "el": atom.GetSymbol(),
                "x": round(p.x, 6),
                "y": round(p.y, 6),
                "z": round(p.z, 6),
                "charge": atom.GetFormalCharge(),
            }
        )
    bonds = []
    odd_bonds = 0
    for bond in work.GetBonds():
        order = _TYPE_TO_ORDER.get(bond.GetBondType())
        if order is None:
            odd_bonds += 1
            order = 1
        bonds.append({"a": bond.GetBeginAtomIdx(), "b": bond.GetEndAtomIdx(), "order": order})
    if odd_bonds and notes is not None:
        notes.append(f"{odd_bonds} bond(s) of an unsupported type were stored as single bonds.")
    return {"atoms": atoms, "bonds": bonds}


def _radical_note(mol: Chem.Mol) -> str | None:
    radicals = [a.GetIdx() for a in mol.GetAtoms() if a.GetNumRadicalElectrons() > 0]
    if not radicals:
        return None
    shown = ", ".join(str(i + 1) for i in radicals[:8]) + (", ..." if len(radicals) > 8 else "")
    return (
        f"{len(radicals)} atom(s) have open valences and were treated as radicals "
        f"(atoms {shown}). Use 'Add hydrogens' first if that is not intended."
    )


# --------------------------------------------------------------------------
# Structure generation and editing helpers
# --------------------------------------------------------------------------


def _embed_3d(mol: Chem.Mol) -> Chem.Mol:
    """Return ``mol`` with explicit hydrogens and a force-field-relaxed 3D conformer."""
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = 0x5EED  # reproducible structures
    if rdDistGeom.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if rdDistGeom.EmbedMolecule(mol, params) != 0:
            raise ChemError("Could not generate 3D coordinates for this structure")
    _relax(mol)
    return mol


def _relax(mol: Chem.Mol) -> str | None:
    """Optimize in place with MMFF94, falling back to UFF. Returns the method used."""
    if FFH.MMFFHasAllMoleculeParams(mol):
        FFH.MMFFOptimizeMolecule(mol, maxIters=2000, ignoreInterfragInteractions=False)
        return "MMFF94"
    if FFH.UFFHasAllMoleculeParams(mol):
        FFH.UFFOptimizeMolecule(mol, maxIters=2000, ignoreInterfragInteractions=False)
        return "UFF"
    return None


def from_smiles(smiles: object) -> tuple[dict, list[str]]:
    """Generate a 3D structure (with hydrogens) from a SMILES string."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ChemError("Enter a SMILES string")
    smiles = smiles.strip()
    if len(smiles) > 2000:
        raise ChemError("SMILES string is too long")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ChemError(f"Could not parse SMILES: {smiles}")
    if mol.GetNumAtoms() == 0:
        raise ChemError("SMILES contains no atoms")
    notes: list[str] = []
    mol = _embed_3d(mol)
    return mol_to_json(mol, notes), notes


def add_hydrogens(data: object) -> tuple[dict, list[str]]:
    """Fill open valences with hydrogens. Existing atoms keep their indices."""
    mol = mol_from_json(data, implicit_h=True)
    before = mol.GetNumAtoms()
    mol = Chem.AddHs(mol, addCoords=True)
    added = mol.GetNumAtoms() - before
    notes = [f"Added {added} hydrogen(s)." if added else "No open valences to fill."]
    return mol_to_json(mol, notes), notes


# --------------------------------------------------------------------------
# Geometry optimization
# --------------------------------------------------------------------------


def optimize(
    data: object,
    method: object = "MMFF94",
    max_iters: object = 2000,
    fixed: object = (),
) -> dict:
    """Minimize the energy of ``data`` with MMFF94, MMFF94s or UFF.

    Returns a dict with the optimized ``molecule``, ``energy`` (kcal/mol, as
    reported by RDKit for the chosen force field), ``initial_energy``,
    ``converged`` and user-facing ``notes``.
    """
    if method not in FORCE_FIELDS:
        raise ChemError(f"Unknown force field {method!r}; choose one of {', '.join(FORCE_FIELDS)}")
    max_iters = _integer(max_iters, "max_iters")
    if not 1 <= max_iters <= 100000:
        raise ChemError("max_iters must be between 1 and 100000")
    mol = mol_from_json(data)
    n = mol.GetNumAtoms()
    if n == 0:
        raise ChemError("Nothing to optimize: the molecule is empty")
    if not isinstance(fixed, (list, tuple)):
        raise ChemError("fixed must be a list of atom indices")
    fixed_atoms = sorted({_integer(i, "Fixed atom index") for i in fixed})
    if any(not 0 <= i < n for i in fixed_atoms):
        raise ChemError("A fixed atom index is out of range")

    notes: list[str] = []
    radical = _radical_note(mol)
    if radical:
        notes.append(radical)

    if method == "UFF":
        if not FFH.UFFHasAllMoleculeParams(mol):
            raise ChemError("UFF has no parameters for some atoms in this molecule.")
        ff = FFH.UFFGetMoleculeForceField(mol, ignoreInterfragInteractions=False)
    else:
        props = FFH.MMFFGetMoleculeProperties(mol, mmffVariant=method)
        if props is None:
            raise ChemError(
                f"{method} has no parameters for this molecule (an element or bonding "
                "situation it does not cover). Try UFF instead."
            )
        ff = FFH.MMFFGetMoleculeForceField(mol, props, ignoreInterfragInteractions=False)
    if ff is None:
        raise ChemError(f"Could not set up the {method} force field for this molecule")

    for i in fixed_atoms:
        ff.AddFixedPoint(i)
    ff.Initialize()
    initial_energy = ff.CalcEnergy()

    # Minimize() restarts the BFGS search on each call; running it in chunks
    # keeps every single call bounded and lets us report convergence.
    converged = False
    done = 0
    while done < max_iters:
        chunk = min(500, max_iters - done)
        if ff.Minimize(maxIts=chunk) == 0:
            converged = True
            break
        done += chunk
    energy = ff.CalcEnergy()

    positions = ff.Positions()
    conf = mol.GetConformer()
    for i in range(n):
        conf.SetAtomPosition(i, Point3D(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]))

    if not converged:
        notes.append(f"Did not fully converge in {max_iters} iterations; run Optimize again to continue.")
    if not all(math.isfinite(c) for c in positions) or not math.isfinite(energy):
        raise ChemError("The optimization diverged (non-finite coordinates). Check the structure.")
    return {
        "molecule": mol_to_json(mol, notes),
        "energy": energy,
        "initial_energy": initial_energy,
        "converged": converged,
        "method": method,
        "notes": notes,
    }


# --------------------------------------------------------------------------
# File import / export
# --------------------------------------------------------------------------


def _extension(name: object) -> str:
    if not isinstance(name, str) or "." not in name:
        raise ChemError("File name needs an extension (.xyz, .mol, .sdf, .mol2 or .pdb)")
    ext = name.rsplit(".", 1)[1].lower()
    if ext not in IMPORT_EXTENSIONS:
        raise ChemError(f"Unsupported file type '.{ext}'. Supported: .xyz, .mol, .sdf, .mol2, .pdb")
    return ext


def import_text(name: object, text: object) -> tuple[dict, list[str]]:
    """Parse a molecule file. Returns the molecule and user-facing notes."""
    ext = _extension(name)
    if not isinstance(text, str) or not text.strip():
        raise ChemError("The file is empty")
    notes: list[str] = []
    if ext == "xyz":
        mol = _read_xyz(text, notes)
    elif ext in ("mol", "sdf", "sd"):
        mol = _read_mol_or_sdf(text, notes)
    elif ext == "mol2":
        mol = Chem.MolFromMol2Block(text, removeHs=False)
        if mol is None:
            raise ChemError("Could not read the Mol2 file")
    else:
        mol = Chem.MolFromPDBBlock(text, removeHs=False)
        if mol is None:
            mol = Chem.MolFromPDBBlock(text, removeHs=False, sanitize=False)
            if mol is None:
                raise ChemError("Could not read the PDB file")
            notes.append("The PDB structure did not pass RDKit's valence checks; loaded as is.")
    if mol.GetNumAtoms() == 0:
        raise ChemError("The file contains no atoms")
    if mol.GetNumConformers() == 0:
        raise ChemError("The file contains no coordinates")
    if ext in ("mol", "sdf", "sd") and not mol.GetConformer().Is3D():
        mol = _embed_3d(mol)
        notes.append("The file had 2D coordinates; generated a 3D structure (hydrogens added).")
    return mol_to_json(mol, notes), notes


def _read_xyz(text: str, notes: list[str]) -> Chem.Mol:
    lines = text.splitlines()
    try:
        count = int(lines[0].split()[0])
    except (IndexError, ValueError):
        raise ChemError("XYZ file must start with the number of atoms") from None
    if count <= 0:
        raise ChemError("XYZ file declares no atoms")
    if len(lines) < count + 2:
        raise ChemError(f"XYZ file declares {count} atoms but is shorter than that")
    frame_lines = lines[: count + 2]
    if any(line.strip() for line in lines[count + 2 :]):
        notes.append("The XYZ file has more than one frame; loaded the first.")
    # Normalize element symbols such as "CL" or "c" and drop extra columns.
    body = [frame_lines[0].strip(), ""]
    for i, line in enumerate(frame_lines[2:]):
        parts = line.split()
        if len(parts) < 4:
            raise ChemError(f"XYZ atom line {i + 1} needs an element and three coordinates")
        sym = parts[0]
        if sym.isdigit():
            z = int(sym)
            if not 0 < z < 119:
                raise ChemError(f"XYZ atom line {i + 1}: invalid atomic number {sym}")
            sym = _PT.GetElementSymbol(z)
        else:
            sym = sym[:1].upper() + sym[1:].lower()
        _atomic_number(sym)
        coords = [_parse_float(p, i + 1) for p in parts[1:4]]
        body.append(f"{sym} {coords[0]:.8f} {coords[1]:.8f} {coords[2]:.8f}")
    mol = Chem.MolFromXYZBlock("\n".join(body) + "\n")
    if mol is None:
        raise ChemError("Could not read the XYZ file")
    if mol.GetNumAtoms() > MAX_ATOMS_FOR_BOND_ORDERS:
        rdDetermineBonds.DetermineConnectivity(mol)
        notes.append(
            f"More than {MAX_ATOMS_FOR_BOND_ORDERS} atoms: connectivity was perceived but all "
            "bonds were set to single."
        )
        return mol
    try:
        work = Chem.Mol(mol)
        rdDetermineBonds.DetermineBonds(work, charge=0)
        return work
    except Exception:
        rdDetermineBonds.DetermineConnectivity(mol)
        notes.append(
            "Bond orders could not be determined for a neutral molecule; all bonds were set "
            "to single. Fix bond orders and charges by hand if needed."
        )
        return mol


def _parse_float(text: str, line_no: int) -> float:
    try:
        value = float(text.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise ChemError(f"XYZ atom line {line_no}: '{text}' is not a number") from None
    if not math.isfinite(value):
        raise ChemError(f"XYZ atom line {line_no}: coordinates must be finite")
    return value


def _read_mol_or_sdf(text: str, notes: list[str]) -> Chem.Mol:
    blocks = [b for b in re.split(r"^\$\$\$\$[^\n]*$", text, flags=re.M) if b.strip()]
    if not blocks:
        raise ChemError("The file contains no molecule")
    if len(blocks) > 1:
        notes.append(f"The file holds {len(blocks)} molecules; loaded the first.")
    # The first record is used untouched: its first line is a (possibly blank) title.
    mol = Chem.MolFromMolBlock(blocks[0], removeHs=False)
    if mol is None:
        raise ChemError("Could not read the MOL/SDF file (it may have an invalid structure)")
    return mol


def export_text(data: object, fmt: object, title: str = "minimol") -> str:
    """Write a molecule as XYZ, MOL, SDF or PDB text."""
    if fmt not in EXPORT_FORMATS:
        raise ChemError(f"Unsupported export format {fmt!r}")
    validate_json(data)
    if fmt == "xyz":
        atoms = data["atoms"]
        lines = [str(len(atoms)), title]
        for a in atoms:
            lines.append(f"{a['el']:<2} {a['x']:14.8f} {a['y']:14.8f} {a['z']:14.8f}")
        return "\n".join(lines) + "\n"
    try:
        mol = mol_from_json(data)
    except ChemError:
        # Still allow saving structures that are chemically odd mid-edit.
        mol = mol_from_json(data, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
    mol.SetProp("_Name", title)
    if fmt == "pdb":
        return Chem.MolToPDBBlock(mol)
    block = Chem.MolToMolBlock(mol, kekulize=False)
    return block + "$$$$\n" if fmt == "sdf" else block
