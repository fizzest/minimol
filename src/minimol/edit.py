"""Structure edits with hydrogen bookkeeping, fragment hotkeys and ring fusion.

Every edit takes a molecule in minimol JSON, applies one operation and returns

    {"molecule": ..., "map": [...], "focus": int | None, "notes": [...]}

where ``map[old_index]`` is the atom's new index (or -1 if it was removed) and
``focus`` is the atom the editor should treat as the new "hotspot" (for
example the group just added), so key sequences can continue from it.

Edits work on 3D molecules (explicit hydrogens) and on 2D drawings
(``dim=2``: flat, uniform bond length, implicit hydrogens, ChemDraw-style
placement; see layout.py).

Hydrogen rules (3D)
-------------------
* Bonding to an atom that has no free valence replaces one of its hydrogens
  (the one pointing most towards where the new bond goes).
* With ``keep_h=True`` the hydrogens of every atom an edit touches are
  re-derived from RDKit's valence model. Hydrogens on untouched atoms are
  never moved. With ``keep_h=False`` nothing is added.

Fragment hotkeys (hover an atom, press a key; modelled on ChemDraw)
------------------------------------------------------------------
* "1" adds a carbon to the hovered atom ("4"/"5": a methyl towards / away from
  the viewer, the 3D meaning of a wedge / hash).
* Every other key turns the hovered atom into the group: "3" on a CH3 makes it
  a phenyl ring, "o" makes it OH. On a hydrogen, that hydrogen becomes the
  group. On an atom inside a chain or ring, single-atom groups change the
  element and larger groups are added as a substituent.

Bond hotkeys: "1"/"2"/"3" set the order, "a" fuses a benzene ring, "v" and
"4"-"8" fuse a 3- to 8-membered ring on the open side of the bond.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdDistGeom

from . import layout
from .chem import ChemError, _atomic_number, _integer, _relax, mol_from_json, mol_to_json

_PT = Chem.GetPeriodicTable()
ORIG = "_mm_orig"
FOCUS = "_mm_focus"
ADJUST = "_mm_adjust"
MAX_CHARGE = 8
FLAT_BOND = layout.BOND_LENGTH


# --------------------------------------------------------------------------
# Fragment library
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Fragment:
    key: str
    name: str
    smiles: str  # one "*" marks the attachment point; atom map 1 marks the focus atom
    mode: str  # "attach": bond a new group to the hovered atom; "replace": the hovered atom becomes the group
    side: str = ""  # "front" / "back": attach towards / away from the viewer


FRAGMENTS = (
    # Keys that add a new atom to the hovered one.
    Fragment("1", "carbon (chain)", "*C", "attach"),
    Fragment("4", "methyl, towards you", "*C", "attach", "front"),
    Fragment("5", "methyl, away from you", "*C", "attach", "back"),
    # Keys that turn the hovered atom into a group (it becomes the attachment atom).
    Fragment("2", "C=O (carbonyl C)", "*C=O", "replace"),
    Fragment("3", "phenyl", "*c1cc[cH:1]cc1", "replace"),
    Fragment("a", "phenyl", "*c1cc[cH:1]cc1", "replace"),
    Fragment("6", "cyclohexyl", "*C1CC[CH2:1]CC1", "replace"),
    Fragment("7", "cyclopentyl", "*C1C[CH2:1]CC1", "replace"),
    Fragment("8", "vinyl", "*C=[CH2:1]", "replace"),
    Fragment("9", "isopropyl", "*C(C)C", "replace"),
    Fragment("z", "ethynyl", "*C#[CH:1]", "replace"),
    Fragment("Z", "azide", "*N=[N+]=[N-]", "replace"),
    Fragment("v", "cyclopropyl", "*C1C[CH2:1]1", "replace"),
    Fragment("u", "cyclobutyl", "*C1C[CH2:1]C1", "replace"),
    Fragment("k", "SO2", "*S(=O)=O", "replace"),
    Fragment("K", "tert-butyl", "*C(C)(C)C", "replace"),
    Fragment("o", "OH", "*O", "replace"),
    Fragment("q", "OH", "*O", "replace"),
    Fragment("O", "OMe", "*OC", "replace"),
    Fragment("n", "NH2", "*N", "replace"),
    Fragment("w", "NH2", "*N", "replace"),
    Fragment("N", "NO2", "*[N+](=O)[O-]", "replace"),
    Fragment("s", "SH", "*S", "replace"),
    Fragment("S", "SiH3", "*[SiH3]", "replace"),
    Fragment("p", "PH2", "*P", "replace"),
    Fragment("P", "Ph", "*c1cc[cH:1]cc1", "replace"),
    Fragment("f", "F", "*F", "replace"),
    Fragment("F", "CF3", "*C(F)(F)F", "replace"),
    Fragment("l", "Cl", "*Cl", "replace"),
    Fragment("C", "Cl", "*Cl", "replace"),
    Fragment("L", "Li", "*[Li]", "replace"),
    Fragment("b", "Br", "*Br", "replace"),
    Fragment("B", "BH2", "*B", "replace"),
    Fragment("i", "I", "*I", "replace"),
    Fragment("h", "H", "*[H]", "replace"),
    Fragment("m", "CH3", "*C", "replace"),
    Fragment("c", "CH3", "*C", "replace"),
    Fragment("e", "Et", "*CC", "replace"),
    Fragment("A", "Ac", "*C(C)=O", "replace"),
    Fragment("E", "CO2Me", "*C(=O)OC", "replace"),
    Fragment("M", "MgBr", "*[Mg]Br", "replace"),
    Fragment("y", "Boc", "*C(=O)OC(C)(C)C", "replace"),
    Fragment("H", "Cbz", "*C(=O)OCc1ccccc1", "replace"),
    Fragment("Q", "Fmoc", "*C(=O)OCC1c2ccccc2-c2ccccc21", "replace"),
)
FRAGMENT_BY_KEY = {f.key: f for f in FRAGMENTS}

# Bond hotkeys: key -> ("order", n) or ("fuse", ring size, aromatic)
BOND_KEYS = {
    "1": ("order", 1),
    "2": ("order", 2),
    "3": ("order", 3),
    "a": ("fuse", 6, True),
    "v": ("fuse", 3, False),
    "4": ("fuse", 4, False),
    "5": ("fuse", 5, False),
    "6": ("fuse", 6, False),
    "7": ("fuse", 7, False),
    "8": ("fuse", 8, False),
}
BOND_KEY_NAMES = {
    "1": "single", "2": "double", "3": "triple", "a": "fuse benzene", "v": "fuse cyclopropane",
    "4": "fuse cyclobutane", "5": "fuse cyclopentane", "6": "fuse cyclohexane",
    "7": "fuse cycloheptane", "8": "fuse cyclooctane",
}


@dataclass(frozen=True)
class Template:
    elements: tuple[str, ...]
    charges: tuple[int, ...]
    coords: np.ndarray  # (n, 3)
    bonds: tuple[tuple[int, int, int], ...]
    dummy: int  # stands in for the parent atom
    attach: int
    focus: int

    @property
    def heavy_count(self) -> int:
        return sum(1 for i, e in enumerate(self.elements) if e != "H" and i != self.dummy)


_ORDER_OF = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}


@lru_cache(maxsize=None)
def template(key: str, flat: bool = False) -> Template:
    """Geometry of a fragment, with an H (3D) or a dummy (2D) in place of the attachment point."""
    frag = FRAGMENT_BY_KEY[key]
    mol = Chem.MolFromSmiles(frag.smiles, sanitize=False)
    dummies = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) != 1:
        raise ValueError(f"fragment {key!r} needs exactly one attachment point")
    dummy = dummies[0]
    nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(dummy).GetNeighbors()]
    if len(nbrs) != 1:
        raise ValueError(f"fragment {key!r}: the attachment point must have one neighbour")
    attach = nbrs[0]
    focus = attach
    rw = Chem.RWMol(mol)
    for a in rw.GetAtoms():
        if a.GetAtomMapNum() == 1:
            focus = a.GetIdx()
        a.SetAtomMapNum(0)
    rw.GetAtomWithIdx(dummy).SetAtomicNum(1)  # a hydrogen stands in for the parent
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    if flat:
        rdDepictor.Compute2DCoords(m)  # bonds of length 1.5, like layout.BOND_LENGTH
    else:
        m = Chem.AddHs(m)
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 0x5EED
        if rdDistGeom.EmbedMolecule(m, params) != 0:
            raise ValueError(f"could not embed fragment {key!r}")
        _relax(m)
    Chem.Kekulize(m, clearAromaticFlags=True)
    conf = m.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    if flat:
        coords[:, 2] = 0.0
    return Template(
        elements=tuple(a.GetSymbol() for a in m.GetAtoms()),
        charges=tuple(a.GetFormalCharge() for a in m.GetAtoms()),
        coords=coords,
        bonds=tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), _ORDER_OF[b.GetBondType()]) for b in m.GetBonds()),
        dummy=dummy,
        attach=attach,
        focus=focus,
    )


def fragment_list() -> list[dict]:
    return [{"key": f.key, "name": f.name, "mode": f.mode} for f in FRAGMENTS]


def bond_key_list() -> list[dict]:
    out = []
    for key, spec in BOND_KEYS.items():
        if spec[0] == "order":
            out.append({"key": key, "name": BOND_KEY_NAMES[key], "op": "order", "order": spec[1]})
        else:
            out.append({"key": key, "name": BOND_KEY_NAMES[key], "op": "fuse", "size": spec[1], "aromatic": spec[2]})
    return out


# --------------------------------------------------------------------------
# Small vector helpers
# --------------------------------------------------------------------------


def _unit(v: np.ndarray) -> np.ndarray | None:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else None


def _align(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector u onto unit vector v."""
    c = float(np.dot(u, v))
    if c > 1 - 1e-12:
        return np.eye(3)
    if c < -1 + 1e-12:
        # Half turn; for vectors in the xy-plane turn about z so they stay in it.
        if abs(u[2]) < 1e-9 and abs(v[2]) < 1e-9:
            return _axis_rotation(np.array([0, 0, 1.0]), math.pi)
        other = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
        return _axis_rotation(_unit(np.cross(u, other)), math.pi)
    k = np.cross(u, v)
    s = float(np.linalg.norm(k))
    k = k / s
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + s * kx + (1 - c) * kx @ kx


def _axis_rotation(k: np.ndarray, angle: float) -> np.ndarray:
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(angle) * kx + (1 - math.cos(angle)) * kx @ kx


def _perpendicular(u: np.ndarray, hint: np.ndarray) -> np.ndarray:
    p = hint - u * float(np.dot(hint, u))
    p = _unit(p)
    if p is None:
        other = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
        p = _unit(np.cross(u, other))
    return p


def _rot2d(v: np.ndarray, degrees: float) -> np.ndarray:
    t = math.radians(degrees)
    return np.array([v[0] * math.cos(t) - v[1] * math.sin(t), v[0] * math.sin(t) + v[1] * math.cos(t), 0.0])


# --------------------------------------------------------------------------
# Editing session on an RDKit molecule
# --------------------------------------------------------------------------


@dataclass
class View:
    right: np.ndarray
    up: np.ndarray
    to_camera: np.ndarray


FLAT_VIEW = View(np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))


def _vector(value: object, what: str) -> np.ndarray:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ChemError(f"{what} must be a list of three numbers")
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in value):
        raise ChemError(f"{what} must contain finite numbers")
    return np.array([float(x) for x in value])


def _view(value: object) -> View:
    if value is None:
        return FLAT_VIEW
    if not isinstance(value, dict):
        raise ChemError("view must be an object")
    right = _unit(_vector(value.get("right"), "view.right"))
    up = _unit(_vector(value.get("up"), "view.up"))
    to_camera = _unit(_vector(value.get("toCamera"), "view.toCamera"))
    if right is None or up is None or to_camera is None:
        raise ChemError("view axes must not be zero")
    return View(right, up, to_camera)


class Session:
    """Mutable working copy used while applying one edit."""

    def __init__(self, data: object, keep_h: bool, view: View, flat: bool = False):
        self.mol = mol_from_json(data, sanitize=False)
        for a in self.mol.GetAtoms():
            a.SetIntProp(ORIG, a.GetIdx())
        self.n_orig = self.mol.GetNumAtoms()
        self.flat = flat
        self.keep_h = keep_h and not flat  # 2D drawings keep hydrogens implicit
        self.view = FLAT_VIEW if flat else view
        self.notes: list[str] = []
        self.doomed: set[int] = set()  # atoms to remove at the end (current indices, nothing removed yet)

    # ---------------------------------------------------------------- queries

    def check_atom(self, i: object, what: str = "atom") -> int:
        i = _integer(i, what)
        if not 0 <= i < self.mol.GetNumAtoms():
            raise ChemError(f"{what} {i + 1} does not exist")
        return i

    def pos(self, i: int) -> np.ndarray:
        return np.array(list(self.mol.GetConformer().GetAtomPosition(i)))

    def el(self, i: int) -> str:
        return self.mol.GetAtomWithIdx(i).GetSymbol()

    def is_h(self, i: int) -> bool:
        return self.mol.GetAtomWithIdx(i).GetAtomicNum() == 1

    def neighbors(self, i: int) -> list[int]:
        return [n.GetIdx() for n in self.mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in self.doomed]

    def heavy_neighbors(self, i: int) -> list[int]:
        return [j for j in self.neighbors(i) if not self.is_h(j)]

    def h_neighbors(self, i: int) -> list[int]:
        """Hydrogens bonded only to atom i."""
        return [j for j in self.neighbors(i) if self.is_h(j) and self.mol.GetAtomWithIdx(j).GetDegree() == 1]

    def open_valence(self, i: int) -> int:
        """How many hydrogens RDKit's valence model would still put on atom i."""
        probe = Chem.RWMol(self.mol)
        for idx in sorted(self.doomed, reverse=True):
            probe.RemoveAtom(idx)
        j = i - sum(1 for d in self.doomed if d < i)
        for a in probe.GetAtoms():
            a.SetNoImplicit(True)
        probe.GetAtomWithIdx(j).SetNoImplicit(False)
        probe.UpdatePropertyCache(strict=False)
        return probe.GetAtomWithIdx(j).GetNumImplicitHs()

    def bond_length(self, a: str, b: str) -> float:
        return FLAT_BOND if self.flat else _PT.GetRcovalent(a) + _PT.GetRcovalent(b)

    def clash(self, p: np.ndarray, ignore: set[int]) -> float:
        """Distance from p to the closest atom not in ``ignore``."""
        best = math.inf
        for j in range(self.mol.GetNumAtoms()):
            if j not in self.doomed and j not in ignore:
                best = min(best, float(np.linalg.norm(self.pos(j) - p)))
        return best

    def default_hint(self, i: int) -> np.ndarray:
        """Away from the heavy neighbours, within the screen plane."""
        p = self.pos(i)
        away = np.zeros(3)
        for j in self.heavy_neighbors(i):
            u = _unit(self.pos(j) - p)
            if u is not None:
                away -= u
        away -= self.view.to_camera * float(np.dot(away, self.view.to_camera))
        return _unit(away) if _unit(away) is not None else self.view.right

    def auto_direction(self, i: int, hint: np.ndarray | None) -> np.ndarray:
        """Direction for a new bond on atom i when no hydrogen is replaced."""
        if self.flat:
            return self.flat_direction(i, hint)
        p = self.pos(i)
        dirs = [u for u in (_unit(self.pos(j) - p) for j in self.neighbors(i)) if u is not None]
        if hint is not None:
            return hint
        if not dirs:
            return self.view.right
        if len(dirs) == 1:
            u = dirs[0]
            w = _perpendicular(u, self.view.up)
            t = math.radians(109.47)
            return _unit(u * math.cos(t) + w * math.sin(t))
        s = np.sum(dirs, axis=0)
        if np.linalg.norm(s) > 0.2:
            return _unit(-s)
        n = np.cross(dirs[0], dirs[1])
        n = _unit(n) if np.linalg.norm(n) > 0.1 else _perpendicular(dirs[0], self.view.up)
        return -n if float(np.dot(n, self.view.to_camera)) < 0 else n

    def flat_direction(self, i: int, hint: np.ndarray | None) -> np.ndarray:
        """ChemDraw-like placement in 2D: 30 degrees for a first bond, 120 degree zigzag, largest gap."""
        if hint is not None:
            h = _unit(np.array([hint[0], hint[1], 0.0]))
            if h is not None:
                return h
        p = self.pos(i)
        nbs = [j for j in self.neighbors(i) if _unit(self.pos(j) - p) is not None]
        if not nbs:
            return _rot2d(np.array([1.0, 0, 0]), 30)
        dirs = [_unit(self.pos(j) - p) for j in nbs]
        if len(dirs) == 1:
            u = dirs[0]
            candidates = [_rot2d(u, 120), _rot2d(u, -120)]
            parent = nbs[0]
            grand = [k for k in self.neighbors(parent) if k != i]
            if grand:
                # Zigzag: the new bond runs parallel to the bond grandparent -> parent.
                want = _unit(self.pos(parent) - self.pos(grand[0]))
                score = [float(np.dot(c, want)) for c in candidates]
            else:
                score = [-abs(c[1]) for c in candidates]  # prefer the more horizontal one
            order = sorted(range(2), key=lambda k: -score[k])
            for k in order:
                if self.clash(p + candidates[k] * FLAT_BOND, {i}) > 0.6 * FLAT_BOND:
                    return candidates[k]
            return candidates[order[0]]
        angles = sorted(math.atan2(d[1], d[0]) for d in dirs)
        gaps = [(angles[(k + 1) % len(angles)] - angles[k]) % (2 * math.pi) or 2 * math.pi for k in range(len(angles))]
        k = max(range(len(gaps)), key=lambda x: gaps[x])
        mid = angles[k] + gaps[k] / 2
        return np.array([math.cos(mid), math.sin(mid), 0.0])

    def slot_for_new_bond(self, i: int, hint: np.ndarray | None) -> np.ndarray:
        """Direction for a new bond on atom i, replacing a hydrogen if i has no free valence.

        With a ``hint`` (a drag direction) the new bond follows the hint exactly,
        so it lands where the editor's preview showed it; the hydrogen closest to
        that direction is the one replaced. Without a hint the new bond takes
        the replaced hydrogen's exact place.
        """
        hs = self.h_neighbors(i)
        if hs and self.open_valence(i) == 0:
            aim = hint if hint is not None else self.default_hint(i)
            p = self.pos(i)
            best = max(hs, key=lambda h: float(np.dot(_unit(self.pos(h) - p), aim)))
            self.doomed.add(best)
            if hint is not None:
                return self.auto_direction(i, hint)
            return _unit(self.pos(best) - p)
        return self.auto_direction(i, hint)

    # ----------------------------------------------------------------- edits

    def add_atom(self, el: str, p: np.ndarray, charge: int = 0) -> int:
        a = Chem.Atom(_atomic_number(el))
        a.SetFormalCharge(charge)
        a.SetNoImplicit(True)
        i = self.mol.AddAtom(a)
        if self.flat:
            p = np.array([p[0], p[1], 0.0])
        self.mol.GetConformer().SetAtomPosition(i, p.tolist())
        return i

    def add_bond(self, a: int, b: int, order: int = 1) -> None:
        self.mol.AddBond(a, b, _ORDER[order])

    def set_focus(self, i: int | None) -> None:
        for a in self.mol.GetAtoms():
            if a.HasProp(FOCUS):
                a.ClearProp(FOCUS)
        if i is not None:
            self.mol.GetAtomWithIdx(i).SetIntProp(FOCUS, 1)

    def mark_adjust(self, i: int) -> None:
        if self.keep_h and i not in self.doomed:
            self.mol.GetAtomWithIdx(i).SetIntProp(ADJUST, 1)

    def doom_with_hydrogens(self, i: int) -> None:
        self.doomed.add(i)
        for h in self.h_neighbors(i):
            self.doomed.add(h)

    def place_template(self, tpl: Template, parent: int | None, base: np.ndarray, d: np.ndarray) -> int:
        """Add a fragment whose attach atom sits at ``base`` with its open bond pointing back along -d.

        The fragment is turned about the new bond (3D) or mirrored (2D) to keep
        clear of existing atoms. Returns the index of the fragment's focus atom.
        """
        A = tpl.coords[tpl.attach]
        u = _unit(tpl.coords[tpl.dummy] - A)
        rot = _align(u, -d)
        local = (tpl.coords - A) @ rot.T  # attach atom at origin, open bond along -d
        keep = [
            k
            for k, e in enumerate(tpl.elements)
            if k != tpl.dummy and (self.keep_h or e != "H" or k == tpl.attach)
        ]
        others = [self.pos(j) for j in range(self.mol.GetNumAtoms()) if j not in self.doomed and j != parent]
        if self.flat:
            # The two orientations in the plane: as is, or mirrored across the bond axis.
            mirror = 2 * np.outer(d, d) - np.eye(3)
            options = [np.eye(3), mirror]
        else:
            options = [_axis_rotation(d, step * math.pi / 12) for step in range(24)]
        best, best_score = options[0], None
        if others and len(keep) > 1:
            existing = np.array(others)
            for option in options:
                pts = local[keep] @ option.T + base
                diff = pts[:, None, :] - existing[None, :, :]
                r2 = np.maximum(np.einsum("ijk,ijk->ij", diff, diff), 0.25)
                score = float(np.sum(1.0 / r2**3))
                if best_score is None or score < best_score - 1e-12:
                    best, best_score = option, score
        final = local @ best.T + base
        index = {}
        for k in keep:
            index[k] = self.add_atom(tpl.elements[k], final[k], tpl.charges[k])
        for a, b, order in tpl.bonds:
            if a in index and b in index:
                self.add_bond(index[a], index[b], order)
        if parent is not None:
            self.add_bond(parent, index[tpl.attach])
        elif self.keep_h:
            # No parent: the open attachment point gets a hydrogen.
            self.mol.GetAtomWithIdx(index[tpl.attach]).SetIntProp(ADJUST, 1)
        return index.get(tpl.focus, index[tpl.attach])

    # --------------------------------------------------------------- finish

    def finish(self) -> dict:
        mol = self.mol
        for idx in sorted(self.doomed, reverse=True):
            mol.RemoveAtom(idx)
        mol = _adjust_hydrogens(mol)
        mapping = [-1] * self.n_orig
        focus = None
        for a in mol.GetAtoms():
            if a.HasProp(ORIG):
                mapping[a.GetIntProp(ORIG)] = a.GetIdx()
            if a.HasProp(FOCUS):
                focus = a.GetIdx()
        if self.flat:
            molecule = layout.to_json_2d(mol, self.notes)  # flags such atoms as "bad" for the drawing
        else:
            self.notes.extend(_valence_notes(mol))
            molecule = mol_to_json(mol, self.notes)
        return {"molecule": molecule, "map": mapping, "focus": focus, "notes": self.notes}


def _valence_notes(mol: Chem.Mol) -> list[str]:
    """Warnings for atoms with more bonds than any allowed valence (3D: hydrogens are explicit)."""
    probe = Chem.RWMol(mol)
    for a in probe.GetAtoms():
        a.SetNoImplicit(True)
    probe.UpdatePropertyCache(strict=False)
    bad = sorted(
        p.GetAtomIdx()
        for p in Chem.DetectChemistryProblems(probe)
        if p.GetType() == "AtomValenceException"
    )
    if not bad:
        return []
    names = ", ".join(f"{probe.GetAtomWithIdx(i).GetSymbol()}{i + 1}" for i in bad[:6]) + (", ..." if len(bad) > 6 else "")
    return [f"Too many bonds on {names}; fix this before optimizing."]


def _adjust_hydrogens(mol: Chem.RWMol) -> Chem.RWMol:
    """Re-derive the hydrogens of atoms flagged with ADJUST (and clear the flag)."""
    flagged = [a.GetIdx() for a in mol.GetAtoms() if a.HasProp(ADJUST)]
    if not flagged:
        return mol
    flagged_set = set(flagged)
    remove = sorted(
        {
            n.GetIdx()
            for i in flagged
            for n in mol.GetAtomWithIdx(i).GetNeighbors()
            if n.GetAtomicNum() == 1 and n.GetDegree() == 1 and n.GetIdx() not in flagged_set
        },
        reverse=True,
    )
    for idx in remove:
        mol.RemoveAtom(idx)
    targets = [a.GetIdx() for a in mol.GetAtoms() if a.HasProp(ADJUST) and a.GetAtomicNum() != 1]
    for a in mol.GetAtoms():
        a.SetNoImplicit(a.GetIdx() not in targets)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    # Without sanitization atoms have no hybridization, and AddHs then puts
    # some of the new hydrogens at the origin. Assign it explicitly.
    Chem.SetHybridization(mol)
    # Hydrogen counts are read before AddHs, which returns a new molecule;
    # existing atoms keep their indices (new hydrogens are appended).
    isolated = {
        i: mol.GetAtomWithIdx(i).GetNumImplicitHs()
        for i in targets
        if mol.GetAtomWithIdx(i).GetDegree() == 0 and mol.GetAtomWithIdx(i).GetNumImplicitHs() > 0
    }
    bonded = [i for i in targets if i not in isolated]
    if bonded:
        mol = Chem.RWMol(Chem.AddHs(mol, addCoords=True, onlyOnAtoms=tuple(bonded)))
    for i, count in isolated.items():
        _cap_isolated(mol, i, count)
    for a in mol.GetAtoms():
        a.SetNoImplicit(True)
        if a.HasProp(ADJUST):
            a.ClearProp(ADJUST)
    return mol


def _cap_isolated(mol: Chem.RWMol, i: int, count: int) -> None:
    """Give an atom without neighbours its hydrogens in a proper geometry.

    RDKit's AddHs(addCoords=True) does not place hydrogens well around an atom
    with no neighbours, so the small molecule (e.g. CH4, NH4+) is embedded
    instead and its hydrogen offsets are copied.
    """
    atom = mol.GetAtomWithIdx(i)
    small = Chem.RWMol()
    a = Chem.Atom(atom.GetAtomicNum())
    a.SetFormalCharge(atom.GetFormalCharge())
    a.SetNoImplicit(True)
    small.AddAtom(a)
    for _ in range(count):
        h = small.AddAtom(Chem.Atom(1))
        small.AddBond(0, h, Chem.BondType.SINGLE)
    small = small.GetMol()
    try:
        Chem.SanitizeMol(small)
        sanitized = True
    except Exception:  # unusual charge/valence: embed without force-field relaxation
        small.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(small)
        sanitized = False
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = 0x5EED
    center = np.array(list(mol.GetConformer().GetAtomPosition(i)))
    if rdDistGeom.EmbedMolecule(small, params) == 0:
        if sanitized:
            _relax(small)
        conf = small.GetConformer()
        origin = np.array(list(conf.GetAtomPosition(0)))
        offsets = [np.array(list(conf.GetAtomPosition(k))) - origin for k in range(1, count + 1)]
    else:  # fall back to a simple spread along the axes
        axes = [np.array(v, float) for v in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))]
        offsets = [axes[k % 6] for k in range(count)]
    length = _PT.GetRcovalent(atom.GetSymbol()) + _PT.GetRcovalent("H")
    for off in offsets:
        h = Chem.Atom(1)
        h.SetNoImplicit(True)
        hi = mol.AddAtom(h)
        mol.AddBond(i, hi, Chem.BondType.SINGLE)
        mol.GetConformer().SetAtomPosition(hi, (center + _unit(off) * length).tolist())


# --------------------------------------------------------------------------
# Public entry point
# --------------------------------------------------------------------------

_ORDER = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}


def apply_edit(data: object, op: object, keep_h: object = True, view: object = None, dim: object = 3) -> dict:
    """Apply one editing operation. See the module docstring for the result."""
    if not isinstance(op, dict):
        raise ChemError("op must be an object")
    if not isinstance(keep_h, bool):
        raise ChemError("keep_h must be true or false")
    if dim not in (2, 3):
        raise ChemError("dim must be 2 or 3")
    s = Session(data, keep_h, _view(view), flat=dim == 2)
    kind = op.get("type")
    hint = _unit(_vector(op["hint"], "hint")) if op.get("hint") is not None else None

    if kind == "add_atom":
        el = op.get("element")
        _atomic_number(el)
        i = s.add_atom(el, _vector(op.get("pos"), "pos"))
        s.mark_adjust(i)
        s.set_focus(i)

    elif kind == "grow":
        i = s.check_atom(op.get("atom"))
        el = op.get("element")
        _atomic_number(el)
        if s.is_h(i) and s.neighbors(i):
            raise ChemError("A hydrogen already has its bond. Click the H with another element to replace it.")
        d = s.slot_for_new_bond(i, hint)
        j = s.add_atom(el, s.pos(i) + d * s.bond_length(s.el(i), el))
        s.add_bond(i, j)
        s.mark_adjust(j)
        s.set_focus(j)

    elif kind == "bond":
        a = s.check_atom(op.get("a"))
        b = s.check_atom(op.get("b"))
        if a == b:
            raise ChemError("An atom cannot bond to itself")
        if s.mol.GetBondBetweenAtoms(a, b) is not None:
            raise ChemError("These atoms are already bonded")
        for x, y in ((a, b), (b, a)):
            if s.is_h(x) and s.neighbors(x):
                raise ChemError("A hydrogen already has its bond")
            s.slot_for_new_bond(x, _unit(s.pos(y) - s.pos(x)))
        s.add_bond(a, b)
        s.set_focus(b)

    elif kind == "element":
        i = s.check_atom(op.get("atom"))
        el = op.get("element")
        s.mol.GetAtomWithIdx(i).SetAtomicNum(_atomic_number(el))
        if el == "H":
            for h in s.h_neighbors(i):
                s.doomed.add(h)
        s.mark_adjust(i)
        s.set_focus(i)

    elif kind == "bond_order":
        a = s.check_atom(op.get("a"))
        b = s.check_atom(op.get("b"))
        order = _integer(op.get("order"), "order")
        bond = s.mol.GetBondBetweenAtoms(a, b)
        if bond is None:
            raise ChemError("These atoms are not bonded")
        if order not in _ORDER:
            raise ChemError("Bond order must be 1, 2 or 3")
        if order > 1 and (s.is_h(a) or s.is_h(b)):
            raise ChemError("Hydrogen can only form a single bond")
        bond.SetBondType(_ORDER[order])
        s.mark_adjust(a)
        s.mark_adjust(b)

    elif kind == "fuse":
        _fuse(s, op, hint)

    elif kind == "charge":
        atoms = op.get("atoms")
        delta = _integer(op.get("delta"), "delta")
        if not isinstance(atoms, list) or not atoms:
            raise ChemError("Choose at least one atom")
        for i in {s.check_atom(x) for x in atoms}:
            atom = s.mol.GetAtomWithIdx(i)
            atom.SetFormalCharge(max(-MAX_CHARGE, min(MAX_CHARGE, atom.GetFormalCharge() + delta)))
            s.mark_adjust(i)

    elif kind == "delete":
        atoms = {s.check_atom(x) for x in op.get("atoms") or []}
        bonds = op.get("bonds") or []
        if not isinstance(bonds, list):
            raise ChemError("bonds must be a list")
        for pair in bonds:
            if not isinstance(pair, list) or len(pair) != 2:
                raise ChemError("Each bond must be a pair of atom indices")
            a, b = (s.check_atom(x) for x in pair)
            if s.mol.GetBondBetweenAtoms(a, b) is None:
                raise ChemError("These atoms are not bonded")
            s.mol.RemoveBond(a, b)
            for x in (a, b):
                if not s.is_h(x):
                    s.mark_adjust(x)
        for i in atoms:
            if s.is_h(i):
                s.doomed.add(i)
            else:
                s.doom_with_hydrogens(i)
        for i in atoms:
            if not s.is_h(i):
                for j in s.heavy_neighbors(i):
                    s.mark_adjust(j)

    elif kind == "fragment":
        _fragment(s, op, hint)

    else:
        raise ChemError(f"Unknown edit {kind!r}")
    return s.finish()


def _fragment(s: Session, op: dict, hint: np.ndarray | None) -> None:
    key = op.get("key")
    frag = FRAGMENT_BY_KEY.get(key) if isinstance(key, str) else None
    if frag is None:
        raise ChemError(f"No fragment for key {key!r}")
    tpl = template(key, s.flat)
    t = s.check_atom(op.get("atom"))
    attach_el = tpl.elements[tpl.attach]

    def place_on(parent: int, d: np.ndarray) -> None:
        base = s.pos(parent) + d * s.bond_length(s.el(parent), attach_el)
        s.set_focus(s.place_template(tpl, parent, base, d))

    def attach_to(parent: int) -> None:
        aim = hint
        if frag.side and not s.flat:
            aim = s.view.to_camera if frag.side == "front" else -s.view.to_camera
        place_on(parent, s.slot_for_new_bond(parent, aim))

    if s.is_h(t):
        nb = s.neighbors(t)
        if len(nb) > 1:
            raise ChemError("This hydrogen has more than one bond")
        s.doomed.add(t)
        if nb:
            place_on(nb[0], _unit(s.pos(t) - s.pos(nb[0])))
        else:
            s.set_focus(s.place_template(tpl, None, s.pos(t), s.view.up))
        return

    if frag.mode == "attach":
        attach_to(t)
        return

    heavy = s.heavy_neighbors(t)
    if len(heavy) == 1:
        parent = heavy[0]
        d = _unit(s.pos(t) - s.pos(parent))
        s.doom_with_hydrogens(t)
        place_on(parent, d)
    elif not heavy:
        p = s.pos(t)
        s.doom_with_hydrogens(t)
        s.set_focus(s.place_template(tpl, None, p, s.view.up))
    elif tpl.heavy_count == 1:
        # An atom inside a chain or ring: change its element instead.
        z = _atomic_number(attach_el)
        valence = sum(
            b.GetBondTypeAsDouble()
            for b in s.mol.GetAtomWithIdx(t).GetBonds()
            if b.GetOtherAtomIdx(t) not in s.doomed and not s.is_h(b.GetOtherAtomIdx(t))
        )
        allowed = max(_PT.GetValenceList(z))
        if allowed != -1 and valence > allowed:
            raise ChemError(f"{attach_el} cannot replace {s.el(t)}{t + 1}: it has {valence:g} bonds to heavy atoms")
        s.mol.GetAtomWithIdx(t).SetAtomicNum(z)
        s.mark_adjust(t)
        s.set_focus(t)
    else:
        # Inside a chain or ring a larger group cannot take the atom's place; add it as a substituent.
        attach_to(t)
        s.notes.append(f"{frag.name} was added to {s.el(t)}{t + 1} as a substituent (the atom is inside a chain or ring).")


def _fuse(s: Session, op: dict, hint: np.ndarray | None) -> None:
    """Fuse a ring onto the bond a-b, like ChemDraw.

    Neighbours of a and b that sit where a ring atom would go (bond angle close
    to the ring's interior angle, with a hydrogen or free valence to give up)
    become part of the new ring: on biphenyl's central bond the ring runs
    through both ortho carbons (phenanthrene); on a ring-CH2 bond it takes in
    the ortho carbon. Atoms of a ring that already contains the bond are never
    reused, so on a benzene edge the ring goes outwards (naphthalene).
    Otherwise all ring atoms are new and the ring goes where there is most room.
    """
    a = s.check_atom(op.get("a"))
    b = s.check_atom(op.get("b"))
    size = _integer(op.get("size"), "size")
    aromatic = op.get("aromatic", False)
    if not 3 <= size <= 8:
        raise ChemError("Ring size must be between 3 and 8")
    if not isinstance(aromatic, bool) or (aromatic and size != 6):
        raise ChemError("Only six-membered rings can be fused as benzene")
    bond = s.mol.GetBondBetweenAtoms(a, b)
    if bond is None:
        raise ChemError("These atoms are not bonded")
    if s.is_h(a) or s.is_h(b):
        raise ChemError("Rings can only be fused onto a bond between two heavy atoms")
    if _ORDER_OF.get(bond.GetBondType(), 1) == 3:
        raise ChemError("A ring cannot be fused onto a triple bond")

    pa_, pb_ = s.pos(a), s.pos(b)
    length = float(np.linalg.norm(pb_ - pa_))
    u = _unit(pb_ - pa_)
    if u is None:
        raise ChemError("The two atoms are on top of each other")
    mid = (pa_ + pb_) / 2
    apothem = length / (2 * math.tan(math.pi / size))
    radius = length / (2 * math.sin(math.pi / size))
    start = math.atan2(-apothem, length / 2)  # angle of b seen from the centre, in the (u, side) basis

    def polygon(side: np.ndarray) -> list[np.ndarray]:
        """Ideal ring vertices after b, going round to the vertex before a (size - 2 points)."""
        center = mid + side * apothem
        return [
            center + radius * (math.cos(start + k * 2 * math.pi / size) * u + math.sin(start + k * 2 * math.pi / size) * side)
            for k in range(1, size - 1)
        ]

    def side_towards(p: np.ndarray) -> np.ndarray | None:
        v = p - mid
        return _unit(v - u * float(np.dot(v, u)))

    partner_a, partner_b = _fusion_partners(s, a, b, size)
    if partner_a is not None or partner_b is not None:
        # Ring through existing neighbours: the side is set by those atoms.
        anchors = [s.pos(x) for x in (partner_a, partner_b) if x is not None]
        side = side_towards(sum(anchors) / len(anchors))
        verts = polygon(side)  # verts[0] is next to b, verts[-1] next to a
        inner = verts[(1 if partner_b is not None else 0) : (len(verts) - 1 if partner_a is not None else len(verts))]
        # Shift the new atoms so they bond sensibly to where the reused atoms really are.
        shift_b = s.pos(partner_b) - verts[0] if partner_b is not None else np.zeros(3)
        shift_a = s.pos(partner_a) - verts[-1] if partner_a is not None else np.zeros(3)
        new = []
        for k, p in enumerate(inner):
            t = (k + 1) / (len(inner) + 1)
            new.append(s.add_atom("C", p + shift_b * (1 - t) + shift_a * t))
        ring = [b, *([partner_b] if partner_b is not None else []), *new, *([partner_a] if partner_a is not None else []), a]
    else:
        new = _fuse_new_atoms(s, a, b, size, aromatic, hint, polygon, u, mid)
        ring = [b, *new, a]

    edges = [(ring[k], ring[(k + 1) % len(ring)]) for k in range(len(ring))]  # closes back to b
    missing = [(x, y) for x, y in edges if s.mol.GetBondBetweenAtoms(x, y) is None]
    if aromatic:
        # New bonds, and existing chain bonds of the new ring, get a Kekule pattern;
        # bonds that already belong to a ring keep their order.
        in_rings = _ring_bonds(s)
        chain = [(x, y) for x, y in edges if (x, y) not in missing and frozenset((x, y)) not in in_rings]
        orders = _kekule_orders(s, ring, missing, chain)
        for (x, y), order in zip(chain, orders[len(missing):]):
            s.mol.GetBondBetweenAtoms(x, y).SetBondType(_ORDER[order])
    else:
        orders = [1] * len(missing)
    for (x, y), order in zip(missing, orders):
        s.add_bond(x, y, order)
    for x in ring:
        s.mark_adjust(x)
    s.set_focus(None)


def _ring_bonds(s: Session) -> set[frozenset]:
    probe = Chem.RWMol(s.mol)
    probe.UpdatePropertyCache(strict=False)
    out = set()
    for ring in Chem.GetSymmSSSR(probe):
        ring = list(ring)
        for k in range(len(ring)):
            out.add(frozenset((ring[k], ring[(k + 1) % len(ring)])))
    return out


def _fusion_partners(s: Session, a: int, b: int, size: int) -> tuple[int | None, int | None]:
    """Neighbours of a and b that the new ring should take in (see _fuse)."""
    interior = 180.0 * (size - 2) / size
    probe = Chem.RWMol(s.mol)
    probe.UpdatePropertyCache(strict=False)
    in_same_ring = set()
    for ring in Chem.GetSymmSSSR(probe):
        if a in ring and b in ring:
            in_same_ring.update(ring)

    def angle(center: int, x: int, y: int) -> float:
        v1, v2 = _unit(s.pos(x) - s.pos(center)), _unit(s.pos(y) - s.pos(center))
        return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2))))))

    def candidates(x: int, other: int) -> list[int]:
        out = []
        for c in s.heavy_neighbors(x):
            if c in (a, b) or c in in_same_ring:
                continue
            if not (s.h_neighbors(c) or s.open_valence(c) > 0):
                continue  # no room for another bond
            if abs(angle(x, other, c) - interior) <= 15:
                out.append(c)
        return out

    cand_a, cand_b = candidates(a, b), candidates(b, a)
    if size >= 5:
        pairs = []
        for pa in cand_a:
            for pb in cand_b:
                if pa == pb or s.mol.GetBondBetweenAtoms(pa, pb) is not None:
                    continue
                dihedral = abs(_dihedral(s.pos(pa), s.pos(a), s.pos(b), s.pos(pb)))
                if dihedral < 70:
                    pairs.append((dihedral, pa, pb))
        if pairs:
            _, pa, pb = min(pairs)
            return pa, pb
    if size >= 4:
        # With one partner, the other end gains a new bond, so it needs room for one
        # (ring - CH2 bond of ethylbenzene: go through the ortho CH, not the CH3).
        def free(x: int) -> bool:
            return bool(s.h_neighbors(x)) or s.open_valence(x) > 0

        singles = [(c, None) for c in cand_a if free(b)] + [(None, c) for c in cand_b if free(a)]
        if singles:
            # Several (e.g. the two ortho carbons): prefer the one with more room around it.
            def room(choice):
                c = choice[0] if choice[0] is not None else choice[1]
                return min(
                    (float(np.linalg.norm(s.pos(c) - s.pos(j))) for j in range(s.mol.GetNumAtoms())
                     if j not in (a, b, c) and not s.is_h(j) and j not in s.heavy_neighbors(c)),
                    default=math.inf,
                )
            return max(singles, key=room)
    return None, None


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    m1 = np.cross(n1, _unit(b1) if _unit(b1) is not None else b1)
    return math.degrees(math.atan2(float(np.dot(m1, n2)), float(np.dot(n1, n2))))


def _fuse_new_atoms(s, a, b, size, aromatic, hint, polygon, u, mid) -> list[int]:
    """All-new ring atoms, on the side of the bond with the most room."""

    def ring_points(side: np.ndarray) -> list[np.ndarray]:
        pts = polygon(side)
        if not s.flat and not aromatic and size >= 5:
            normal = np.cross(u, side)
            # Pucker, so an optimization does not sit on a flat ring.
            pts = [p + normal * (0.25 if k % 2 == 0 else -0.25) for k, p in enumerate(pts)]
        return pts

    base = _perpendicular(u, s.view.up if not s.flat else np.array([-u[1], u[0], 0.0]))
    sides = [base, -base] if s.flat else [_axis_rotation(u, k * math.pi / 12) @ base for k in range(24)]
    # Hydrogens on a and b are re-derived (3D) and do not count as obstacles.
    ignore = {a, b} | (set(s.h_neighbors(a)) | set(s.h_neighbors(b)) if s.keep_h else set())
    others = [s.pos(j) for j in range(s.mol.GetNumAtoms()) if j not in ignore and j not in s.doomed]
    crowd = np.zeros(3)
    for x in (a, b):
        for k in s.heavy_neighbors(x):
            if k not in (a, b):
                v = s.pos(k) - mid
                w = _unit(v - u * float(np.dot(v, u)))
                if w is not None:
                    crowd += w

    def score(side: np.ndarray) -> tuple[float, float]:
        clearance = min((float(np.linalg.norm(p - q)) for p in ring_points(side) for q in others), default=math.inf)
        # Round clearance so near-ties fall back to "away from the neighbours" (then the hint).
        away = -float(np.dot(side, crowd))
        if hint is not None:
            away += 0.01 * float(np.dot(side, hint))
        return (round(min(clearance, 3.0), 1), away)

    side = max(sides, key=score)
    return [s.add_atom("C", p) for p in ring_points(side)]


def _kekule_orders(s: Session, ring: list[int], missing: list[tuple[int, int]], chain: list[tuple[int, int]] = ()) -> list[int]:
    """Orders for the new bonds (``missing``) and the re-assignable ``chain`` bonds of a fused benzene ring.

    Tries every single/double choice and keeps the one that gives the most ring
    atoms exactly one double bond, without two double bonds on an atom or any
    atom over its valence (hydrogens are re-derived afterwards); among equally
    good choices, the one that changes the fewest existing bonds.
    """
    import itertools

    variable = list(missing) + list(chain)
    variable_keys = {frozenset(e) for e in variable}

    def fixed_doubles_and_valence(x: int) -> tuple[int, float]:
        d, v = 0, 0.0
        for bd in s.mol.GetAtomWithIdx(x).GetBonds():
            other = bd.GetOtherAtomIdx(x)
            if s.is_h(other) or frozenset((x, other)) in variable_keys:
                continue
            order = _ORDER_OF.get(bd.GetBondType(), 1)
            d += order == 2
            v += order
        return d, v

    base = {x: fixed_doubles_and_valence(x) for x in ring}
    current = [None] * len(missing) + [_ORDER_OF.get(s.mol.GetBondBetweenAtoms(x, y).GetBondType(), 1) for x, y in chain]
    limit = {}
    for x in ring:
        valences = _PT.GetValenceList(s.mol.GetAtomWithIdx(x).GetAtomicNum())
        limit[x] = math.inf if -1 in valences else max(valences) + s.mol.GetAtomWithIdx(x).GetFormalCharge()
    best, best_score = [1] * len(variable), None
    for choice in itertools.product((2, 1), repeat=len(variable)):
        doubles = {x: base[x][0] for x in ring}
        valence = {x: base[x][1] for x in ring}
        for (x, y), order in zip(variable, choice):
            for z in (x, y):
                doubles[z] += order == 2
                valence[z] += order
        if any(doubles[x] > 1 or valence[x] > limit[x] for x in ring):
            continue
        changed = sum(1 for c, o in zip(current, choice) if c is not None and c != o)
        score = (sum(doubles[x] == 1 for x in ring), -changed)
        if best_score is None or score > best_score:
            best, best_score = list(choice), score
    return best
