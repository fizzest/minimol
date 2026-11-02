"""2D (skeletal, ChemDraw-style) structures and conversion between 2D and 3D.

A 2D molecule uses the same JSON as a 3D one, with every z = 0, bonds of a
uniform length (1.5, RDKit's depiction length) and hydrogens implicit, as in
ChemDraw. Each atom also carries ``h`` (its implicit hydrogen count, for
labels such as "OH") and ``bad`` (true when it exceeds every allowed valence).
These two fields are derived data; the server recomputes them on every 2D
result.
"""

from __future__ import annotations

import statistics

from rdkit import Chem
from rdkit.Chem import rdDepictor

from .chem import ChemError, _embed_3d, mol_from_json, mol_to_json

BOND_LENGTH = 1.5  # RDKit's default depiction bond length


def strip_hydrogens(mol: Chem.Mol) -> Chem.RWMol:
    """Remove hydrogens that hang on one heavy atom; keep H2, bridging or isolated H."""
    rw = Chem.RWMol(mol)
    doomed = [
        a.GetIdx()
        for a in rw.GetAtoms()
        if a.GetAtomicNum() == 1
        and a.GetDegree() == 1
        and a.GetNeighbors()[0].GetAtomicNum() != 1
        and a.GetIsotope() == 0
        and a.GetFormalCharge() == 0
    ]
    for idx in sorted(doomed, reverse=True):
        rw.RemoveAtom(idx)
    return rw


def _prepare(mol: Chem.Mol) -> None:
    for a in mol.GetAtoms():
        a.SetNoImplicit(a.GetAtomicNum() == 1)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)


def depict(mol: Chem.RWMol) -> None:
    """Give ``mol`` fresh 2D coordinates (in place)."""
    _prepare(mol)
    rdDepictor.Compute2DCoords(mol)


def has_2d_coords(mol: Chem.Mol) -> bool:
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    if conf.Is3D() and any(abs(conf.GetAtomPosition(i).z) > 1e-3 for i in range(mol.GetNumAtoms())):
        return False
    return mol.GetNumAtoms() == 1 or any(
        (conf.GetAtomPosition(i) - conf.GetAtomPosition(0)).Length() > 1e-3 for i in range(mol.GetNumAtoms())
    )


def normalize_scale(mol: Chem.Mol) -> None:
    """Scale 2D coordinates so the median bond is BOND_LENGTH long, flatten z (in place)."""
    conf = mol.GetConformer()
    lengths = [
        (conf.GetAtomPosition(b.GetBeginAtomIdx()) - conf.GetAtomPosition(b.GetEndAtomIdx())).Length()
        for b in mol.GetBonds()
    ]
    lengths = [x for x in lengths if x > 1e-6]
    factor = BOND_LENGTH / statistics.median(lengths) if lengths else 1.0
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        p.x *= factor
        p.y *= factor
        p.z = 0.0
        conf.SetAtomPosition(i, p)
    conf.Set3D(False)


def annotate(mol: Chem.Mol, data: dict) -> dict:
    """Add implicit hydrogen counts (``h``) and valence problems (``bad``) to 2D JSON."""
    probe = Chem.RWMol(mol)
    _prepare(probe)
    bad = set()
    for problem in Chem.DetectChemistryProblems(probe):
        if hasattr(problem, "GetAtomIdx") and problem.GetType() == "AtomValenceException":
            bad.add(problem.GetAtomIdx())
    for i, atom in enumerate(data["atoms"]):
        a = probe.GetAtomWithIdx(i)
        atom["h"] = 0 if a.GetAtomicNum() == 1 else a.GetNumImplicitHs()
        atom["z"] = 0.0
        if i in bad:
            atom["bad"] = True
        else:
            atom.pop("bad", None)
    return data


def to_json_2d(mol: Chem.Mol, notes: list[str] | None = None) -> dict:
    return annotate(mol, mol_to_json(mol, notes))


def from_rdkit_2d(mol: Chem.Mol, notes: list[str]) -> dict:
    """A 2D molecule from any RDKit molecule (keeps usable 2D coordinates, else lays it out)."""
    rw = strip_hydrogens(mol)
    if rw.GetNumAtoms() == 0:
        raise ChemError("The structure has no heavy atoms")
    if has_2d_coords(rw):
        normalize_scale(rw)
    else:
        depict(rw)
    return to_json_2d(rw, notes)


def to_2d(data: object) -> dict:
    """Convert a (3D) molecule to a 2D drawing: hydrogens become implicit, fresh layout."""
    mol = strip_hydrogens(mol_from_json(data, sanitize=False))
    if mol.GetNumAtoms() == 0:
        return {"molecule": {"atoms": [], "bonds": []}, "notes": []}
    depict(mol)
    notes: list[str] = []
    return {"molecule": to_json_2d(mol, notes), "notes": notes}


def clean_2d(data: object) -> dict:
    """Lay out a 2D drawing again (ChemDraw's "Clean Structure"), keeping it in place."""
    mol = mol_from_json(data, sanitize=False)
    if mol.GetNumAtoms() == 0:
        return {"molecule": {"atoms": [], "bonds": []}, "notes": []}
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    cx = sum(conf.GetAtomPosition(i).x for i in range(n)) / n
    cy = sum(conf.GetAtomPosition(i).y for i in range(n)) / n
    depict(mol)
    conf = mol.GetConformer()
    nx = sum(conf.GetAtomPosition(i).x for i in range(n)) / n
    ny = sum(conf.GetAtomPosition(i).y for i in range(n)) / n
    for i in range(n):
        p = conf.GetAtomPosition(i)
        p.x += cx - nx
        p.y += cy - ny
        conf.SetAtomPosition(i, p)
    return {"molecule": to_json_2d(mol), "notes": []}


def to_3d(data: object) -> dict:
    """Build a 3D structure (with hydrogens) from a 2D drawing: ETKDG embedding, then MMFF94/UFF.

    Heavy atoms keep their indices; hydrogens are appended.
    """
    mol = mol_from_json(data, sanitize=False)
    if mol.GetNumAtoms() == 0:
        return {"molecule": {"atoms": [], "bonds": []}, "notes": []}
    for a in mol.GetAtoms():
        a.SetNoImplicit(a.GetAtomicNum() == 1)
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        text = str(exc).strip()
        raise ChemError(f"Fix the drawing before switching to 3D: {text}") from None
    mol.RemoveAllConformers()
    notes: list[str] = []
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    unassigned = [a.GetIdx() + 1 for a in mol.GetAtoms() if a.HasProp("_ChiralityPossible")]
    mol3d = _embed_3d(mol)
    if unassigned:
        notes.append(
            "The drawing has no stereo information; stereocentres at atoms "
            + ", ".join(map(str, unassigned[:8]))
            + (", ..." if len(unassigned) > 8 else "")
            + " were set arbitrarily."
        )
    return {"molecule": mol_to_json(mol3d, notes), "notes": notes}


def smiles_2d(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ChemError(f"Could not parse SMILES: {smiles}")
    notes: list[str] = []
    return {"molecule": from_rdkit_2d(mol, notes), "notes": notes}


__all__ = [
    "BOND_LENGTH",
    "annotate",
    "clean_2d",
    "depict",
    "from_rdkit_2d",
    "smiles_2d",
    "strip_hydrogens",
    "to_2d",
    "to_3d",
    "to_json_2d",
]
