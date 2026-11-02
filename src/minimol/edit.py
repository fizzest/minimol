"""Structure edits with hydrogen bookkeeping, and fragment (hotkey) attachment.

Every edit takes a molecule in minimol JSON, applies one operation and returns

    {"molecule": ..., "map": [...], "focus": int | None, "notes": [...]}

where ``map[old_index]`` is the atom's new index (or -1 if it was removed) and
``focus`` is the atom the editor should treat as the new "hotspot" (for
example the atom just added), so key sequences can continue from it.

Hydrogen rules
--------------
* Bonding to an atom that has no free valence replaces one of its hydrogens
  (the one pointing most towards where the new bond goes).
* With ``keep_h=True`` the hydrogens of every atom an edit touches are
  re-derived from RDKit's valence model: new atoms come with hydrogens,
  changing an element, a bond order or a charge fixes that atom's hydrogen
  count, and deleting an atom caps its neighbours with hydrogens. Hydrogens on
  untouched atoms are never moved.
* With ``keep_h=False`` nothing is added; the structure is exactly what is
  drawn (open valences stay open).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom

from .chem import ChemError, _atomic_number, _integer, _relax, mol_from_json, mol_to_json

_PT = Chem.GetPeriodicTable()
ORIG = "_mm_orig"
FOCUS = "_mm_focus"
ADJUST = "_mm_adjust"
MAX_CHARGE = 8


# --------------------------------------------------------------------------
# Fragment library (hover an atom, press the key)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Fragment:
    key: str
    name: str
    smiles: str  # one "*" marks the attachment point; atom map 1 marks the focus atom
    mode: str  # "attach": bond to the hovered atom; "replace": the hovered atom becomes this group


FRAGMENTS = (
    # Build keys: add a group to the hovered atom.
    Fragment("1", "carbon (chain)", "*C", "attach"),
    Fragment("2", "carbonyl C=O", "*C=O", "attach"),
    Fragment("3", "phenyl", "*c1ccccc1", "attach"),
    Fragment("a", "phenyl", "*c1ccccc1", "attach"),
    Fragment("6", "cyclohexyl", "*C1CCCCC1", "attach"),
    Fragment("7", "cyclopentyl", "*C1CCCC1", "attach"),
    Fragment("8", "vinyl", "*C=[CH2:1]", "attach"),
    Fragment("9", "isopropyl", "*C(C)C", "attach"),
    Fragment("z", "ethynyl", "*C#[CH:1]", "attach"),
    Fragment("Z", "azide", "*N=[N+]=[N-]", "attach"),
    Fragment("v", "cyclopropyl", "*C1CC1", "attach"),
    Fragment("u", "cyclobutyl", "*C1CCC1", "attach"),
    Fragment("K", "tert-butyl", "*C(C)(C)C", "attach"),
    # Label keys: the hovered atom (with its hydrogens) becomes this group.
    Fragment("o", "OH", "*O", "replace"),
    Fragment("q", "OH", "*O", "replace"),
    Fragment("O", "OMe", "*OC", "replace"),
    Fragment("n", "NH2", "*N", "replace"),
    Fragment("w", "NH2", "*N", "replace"),
    Fragment("N", "NO2", "*[N+](=O)[O-]", "replace"),
    Fragment("s", "SH", "*S", "replace"),
    Fragment("S", "SiH3", "*[SiH3]", "replace"),
    Fragment("p", "PH2", "*P", "replace"),
    Fragment("P", "Ph", "*c1ccccc1", "replace"),
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


@lru_cache(maxsize=None)
def template(key: str) -> Template:
    """3D geometry of a fragment, with an H in place of the attachment point."""
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
    m = Chem.AddHs(m)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = 0x5EED
    if rdDistGeom.EmbedMolecule(m, params) != 0:
        raise ValueError(f"could not embed fragment {key!r}")
    _relax(m)
    Chem.Kekulize(m, clearAromaticFlags=True)
    conf = m.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    order = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
    return Template(
        elements=tuple(a.GetSymbol() for a in m.GetAtoms()),
        charges=tuple(a.GetFormalCharge() for a in m.GetAtoms()),
        coords=coords,
        bonds=tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), order[b.GetBondType()]) for b in m.GetBonds()),
        dummy=dummy,
        attach=attach,
        focus=focus,
    )


def fragment_list() -> list[dict]:
    return [{"key": f.key, "name": f.name, "mode": f.mode} for f in FRAGMENTS]


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
        other = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
        k = _unit(np.cross(u, other))
        return _axis_rotation(k, math.pi)
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


# --------------------------------------------------------------------------
# Editing session on an RDKit molecule
# --------------------------------------------------------------------------


@dataclass
class View:
    right: np.ndarray
    up: np.ndarray
    to_camera: np.ndarray


def _vector(value: object, what: str) -> np.ndarray:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ChemError(f"{what} must be a list of three numbers")
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in value):
        raise ChemError(f"{what} must contain finite numbers")
    return np.array([float(x) for x in value])


def _view(value: object) -> View:
    if value is None:
        return View(np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
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

    def __init__(self, data: object, keep_h: bool, view: View):
        self.mol = mol_from_json(data, sanitize=False)
        for a in self.mol.GetAtoms():
            a.SetIntProp(ORIG, a.GetIdx())
        self.n_orig = self.mol.GetNumAtoms()
        self.keep_h = keep_h
        self.view = view
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
        return _PT.GetRcovalent(a) + _PT.GetRcovalent(b)

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
            return hint if hint is not None else _unit(self.pos(best) - p)
        if not hs and self.open_valence(i) == 0 and self.heavy_neighbors(i):
            self.notes.append(f"{self.el(i)}{i + 1} had no free valence; it now has an extra bond.")
        return self.auto_direction(i, hint)

    # ----------------------------------------------------------------- edits

    def add_atom(self, el: str, p: np.ndarray, charge: int = 0) -> int:
        a = Chem.Atom(_atomic_number(el))
        a.SetFormalCharge(charge)
        a.SetNoImplicit(True)
        i = self.mol.AddAtom(a)
        self.mol.GetConformer().SetAtomPosition(i, p.tolist())
        return i

    def set_focus(self, i: int) -> None:
        for a in self.mol.GetAtoms():
            if a.HasProp(FOCUS):
                a.ClearProp(FOCUS)
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

        Returns the index of the fragment's focus atom.
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
        others = [
            self.pos(j)
            for j in range(self.mol.GetNumAtoms())
            if j not in self.doomed and j != parent
        ]
        best_angle, best_score = 0.0, None
        if others and len(keep) > 1:
            existing = np.array(others)
            for step in range(24):
                angle = step * math.pi / 12
                pts = local[keep] @ _axis_rotation(d, angle).T + base
                diff = pts[:, None, :] - existing[None, :, :]
                r2 = np.maximum(np.einsum("ijk,ijk->ij", diff, diff), 0.25)
                score = float(np.sum(1.0 / r2**3))
                if best_score is None or score < best_score - 1e-12:
                    best_angle, best_score = angle, score
        final = local @ _axis_rotation(d, best_angle).T + base
        index = {}
        for k in keep:
            index[k] = self.add_atom(tpl.elements[k], final[k], tpl.charges[k])
        for a, b, order in tpl.bonds:
            if a in index and b in index:
                self.mol.AddBond(index[a], index[b], {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[order])
        if parent is not None:
            self.mol.AddBond(parent, index[tpl.attach], Chem.BondType.SINGLE)
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
        return {"molecule": mol_to_json(mol, self.notes), "map": mapping, "focus": focus, "notes": self.notes}


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


def apply_edit(data: object, op: object, keep_h: object = True, view: object = None) -> dict:
    """Apply one editing operation. See the module docstring for the result."""
    if not isinstance(op, dict):
        raise ChemError("op must be an object")
    if not isinstance(keep_h, bool):
        raise ChemError("keep_h must be true or false")
    s = Session(data, keep_h, _view(view))
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
        s.mol.AddBond(i, j, Chem.BondType.SINGLE)
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
        s.mol.AddBond(a, b, Chem.BondType.SINGLE)
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
        _fragment(s, op)

    else:
        raise ChemError(f"Unknown edit {kind!r}")
    return s.finish()


def _fragment(s: Session, op: dict) -> None:
    key = op.get("key")
    frag = FRAGMENT_BY_KEY.get(key) if isinstance(key, str) else None
    if frag is None:
        raise ChemError(f"No fragment for key {key!r}")
    tpl = template(key)
    t = s.check_atom(op.get("atom"))
    hint = _unit(_vector(op["hint"], "hint")) if op.get("hint") is not None else None
    attach_el = tpl.elements[tpl.attach]
    view = s.view

    def place_on(parent: int, d: np.ndarray) -> None:
        base = s.pos(parent) + d * s.bond_length(s.el(parent), attach_el)
        s.set_focus(s.place_template(tpl, parent, base, d))

    if s.is_h(t):
        nb = s.neighbors(t)
        if len(nb) == 1:
            parent = nb[0]
            d = _unit(s.pos(t) - s.pos(parent))
            s.doomed.add(t)
            place_on(parent, d)
            return
        if nb:
            raise ChemError("This hydrogen has more than one bond")
        s.doomed.add(t)
        s.set_focus(s.place_template(tpl, None, s.pos(t), view.up))
        return

    if frag.mode == "attach":
        d = s.slot_for_new_bond(t, hint)
        place_on(t, d)
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
        s.set_focus(s.place_template(tpl, None, p, view.up))
    elif tpl.heavy_count == 1:
        # An atom inside a chain or ring: change its element instead.
        z = _atomic_number(attach_el)
        valence = sum(b.GetBondTypeAsDouble() for b in s.mol.GetAtomWithIdx(t).GetBonds() if b.GetOtherAtomIdx(t) not in s.doomed and not s.is_h(b.GetOtherAtomIdx(t)))
        allowed = max(_PT.GetValenceList(z))
        if allowed != -1 and valence > allowed:
            raise ChemError(f"{attach_el} cannot replace {s.el(t)}{t + 1}: it has {valence:g} bonds to heavy atoms")
        s.mol.GetAtomWithIdx(t).SetAtomicNum(z)
        s.mark_adjust(t)
        s.set_focus(t)
    else:
        raise ChemError(f"{frag.name} replaces end atoms only; hover a terminal atom or a hydrogen")
