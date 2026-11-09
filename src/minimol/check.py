"""Checks of a molecule before a quantum-chemistry run (g-xTB).

g-xTB only sees the atoms' elements, positions and the charge; it knows
nothing about bonds. Three kinds of problem make runs fail or answer a
different question than intended:

* Overlapping or very close atoms: the optimization starts in a nonsense region
  and often does not converge.
* Badly stretched or compressed bonds (for example a template with every bond
  at 1.52 A, also inside aromatic rings): usually fine, but worth knowing.
* An odd number of electrons: no closed-shell state exists, so g-xTB runs an
  open-shell calculation (a doublet) unless told otherwise. Often this is a
  missing or extra hydrogen or a wrong total charge.

The thresholds below are minimol's own rules of thumb, not published limits.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem

from .chem import _atomic_number
from .edit import _valence_for

_PT = Chem.GetPeriodicTable()

# Non-bonded atoms closer than this fraction of the sum of their van der Waals radii are a clash.
CLASH_FRACTION = 0.6
# Non-bonded atoms closer than this (Angstrom) are on top of each other: no run can make sense of that.
OVERLAP = 0.5
# A bond outside this range of the sum of the covalent radii is called stretched / compressed.
BOND_MIN, BOND_MAX = 0.70, 1.30
MAX_LISTED = 4
MAX_ATOMS_FOR_PAIRS = 6000  # the pair check is quadratic


def _name(atoms: list[dict], i: int) -> str:
    return f"{atoms[i]['el']}{i + 1}"


def _list(items: list[str], limit: int = MAX_LISTED) -> str:
    return ", ".join(items[:limit]) + (f" and {len(items) - limit} more" if len(items) > limit else "")


def preflight(data: dict, charge: int, uhf: int | None = None) -> dict:
    """Problems found in ``data`` (a minimol molecule) for a run with total ``charge``.

    Returns ``{"errors": [...], "warnings": [...]}``; atoms are numbered from 1
    like in the editor.
    """
    atoms = data["atoms"]
    bonds = data.get("bonds", [])
    n = len(atoms)
    z = [_atomic_number(a["el"]) for a in atoms]
    errors: list[str] = []
    warnings: list[str] = []
    neighbors: list[set[int]] = [set() for _ in range(n)]
    for b in bonds:
        neighbors[b["a"]].add(b["b"])
        neighbors[b["b"]].add(b["a"])
    xyz = np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float).reshape(n, 3)

    # ---- atoms on top of each other, and clashes
    if 1 < n <= MAX_ATOMS_FOR_PAIRS:
        vdw = np.array([_PT.GetRvdw(x) or 2.0 for x in z])
        overlaps: list[tuple[float, int, int]] = []
        clashes: list[tuple[float, int, int, float]] = []
        for i in range(n - 1):
            d = np.linalg.norm(xyz[i + 1 :] - xyz[i], axis=1)
            limit = CLASH_FRACTION * (vdw[i] + vdw[i + 1 :])
            for k in np.nonzero(d < limit)[0]:
                j = i + 1 + int(k)
                if j in neighbors[i] or neighbors[i] & neighbors[j]:
                    continue  # bonded, or the two ends of an angle: expected to be close
                if d[k] < OVERLAP:
                    overlaps.append((float(d[k]), i, j))
                else:
                    clashes.append((float(d[k]), i, j, float(vdw[i] + vdw[j])))
        if overlaps:
            overlaps.sort()
            errors.append(
                "Atoms on top of each other: "
                + _list([f"{_name(atoms, i)}–{_name(atoms, j)} ({d:.2f} Å)" for d, i, j in overlaps])
                + ". Move them apart (or delete the duplicate) first."
            )
        if clashes:
            clashes.sort()
            warnings.append(
                f"{len(clashes)} pair(s) of non-bonded atoms are far too close, e.g. "
                + _list([f"{_name(atoms, i)}–{_name(atoms, j)} {d:.2f} Å" for d, i, j, _ in clashes], 3)
                + ". g-xTB may not converge from here; run UFF or MMFF94 first to push them apart."
            )

    # ---- bond lengths
    odd_bonds: list[tuple[float, str]] = []
    for b in bonds:
        i, j = b["a"], b["b"]
        expected = _PT.GetRcovalent(z[i]) + _PT.GetRcovalent(z[j])
        length = float(np.linalg.norm(xyz[i] - xyz[j]))
        ratio = length / expected if expected else 1.0
        if ratio > BOND_MAX or ratio < BOND_MIN:
            odd_bonds.append((abs(ratio - 1.0), f"{_name(atoms, i)}–{_name(atoms, j)} {length:.2f} Å (typical {expected:.2f})"))
    if odd_bonds:
        odd_bonds.sort(reverse=True)
        warnings.append(
            f"{len(odd_bonds)} bond(s) have an unusual length, e.g. "
            + _list([text for _, text in odd_bonds], 3)
            + ". Clean the geometry with UFF or MMFF94 before g-xTB."
        )

    # ---- electrons
    electrons = sum(z) - charge
    if electrons % 2 and uhf is None:
        open_valence = []
        for i in range(n):
            if z[i] == 1:
                continue
            used = sum(b.get("order", 1) for b in bonds if i in (b["a"], b["b"]))
            valence = _valence_for(z[i], atoms[i].get("charge", 0), used)
            if valence is not None and valence > used:
                open_valence.append(_name(atoms, i))
        text = (
            f"The electron count is odd ({electrons}), so there is no closed-shell state: "
            "g-xTB will run an open-shell doublet. For a closed-shell molecule change the total charge by ±1 "
            "(+/- buttons) or fix a missing or extra hydrogen"
        )
        if open_valence:
            text += f"; atoms with open valences: {_list(open_valence, 6)}"
        warnings.append(text + ". To keep it open-shell on purpose, enter the unpaired electrons.")
    return {"errors": errors, "warnings": warnings}
