import math

import numpy as np
import pytest
from rdkit import Chem

from minimol import chem, layout
from minimol.chem import ChemError
from minimol.edit import FRAGMENTS, apply_edit


def edit2d(mol, **op):
    return apply_edit(mol, op, keep_h=True, dim=2)


def P(mol, i):
    a = mol["atoms"][i]
    return np.array([a["x"], a["y"], a["z"]])


def bond_lengths(mol):
    return [float(np.linalg.norm(P(mol, b["a"]) - P(mol, b["b"]))) for b in mol["bonds"]]


def flat(mol):
    return all(a["z"] == 0 for a in mol["atoms"])


def smiles_of(mol2d):
    m = chem.mol_from_json(mol2d, sanitize=False)
    for a in m.GetAtoms():
        a.SetNoImplicit(False)
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


EMPTY = {"atoms": [], "bonds": []}


def test_smiles_to_2d_drawing():
    r = layout.smiles_2d("CC(=O)O")
    m = r["molecule"]
    assert len(m["atoms"]) == 4 and flat(m)
    assert all(abs(x - 1.5) < 0.05 for x in bond_lengths(m))
    assert [a["h"] for a in m["atoms"]] == [3, 0, 0, 1]


def test_3d_to_2d_and_back_keeps_heavy_atom_order():
    ethanol, _ = chem.from_smiles("CCO")
    two = layout.to_2d(ethanol)["molecule"]
    assert [a["el"] for a in two["atoms"]] == ["C", "C", "O"]
    assert [a["h"] for a in two["atoms"]] == [3, 2, 1]
    three = layout.to_3d(two)["molecule"]
    assert [a["el"] for a in three["atoms"]][:3] == ["C", "C", "O"]
    assert len(three["atoms"]) == 9 and not flat(three)


def test_to_3d_reports_arbitrary_stereo_and_bad_drawings():
    r = layout.to_3d(layout.smiles_2d("CC(N)C(=O)O")["molecule"])
    assert any("stereo" in n for n in r["notes"])
    bad = {"atoms": [{"el": "C", "x": 0, "y": 0, "z": 0, "charge": 0}] + [{"el": "C", "x": math.cos(k), "y": math.sin(k), "z": 0, "charge": 0} for k in range(5)],
           "bonds": [{"a": 0, "b": k, "order": 1} for k in range(1, 6)]}
    with pytest.raises(ChemError, match="Fix the drawing"):
        layout.to_3d(bad)


def test_2d_drawing_flags_over_valent_atoms():
    m = layout.smiles_2d("CC")["molecule"]
    m = edit2d(m, type="bond_order", a=0, b=1, order=3)["molecule"]  # fine: HC#CH
    assert not any(a.get("bad") for a in m["atoms"])
    m = edit2d(m, type="grow", atom=0, element="C")["molecule"]
    m = edit2d(m, type="grow", atom=0, element="C")["molecule"]  # C with a triple and two singles
    assert m["atoms"][0].get("bad") is True


def test_2d_growth_is_chemdraw_like():
    m = edit2d(EMPTY, type="add_atom", element="C", pos=[0, 0, 0])["molecule"]
    assert len(m["atoms"]) == 1 and m["atoms"][0]["h"] == 4  # no explicit H in 2D
    r = edit2d(m, type="grow", atom=0, element="C")
    m = r["molecule"]
    first = P(m, 1) - P(m, 0)
    assert math.degrees(math.atan2(first[1], first[0])) == pytest.approx(30, abs=1e-4)  # coordinates have 6 decimals
    focus = r["focus"]
    for _ in range(3):
        r = edit2d(m, type="grow", atom=focus, element="C")
        m, focus = r["molecule"], r["focus"]
    assert flat(m) and all(abs(x - 1.5) < 1e-5 for x in bond_lengths(m))
    # Zigzag: every bond angle is 120 degrees and the chain alternates up and down.
    ys = [round(a["y"], 3) for a in m["atoms"]]
    assert len(set(ys)) == 2
    for k in range(1, 4):
        v1, v2 = P(m, k - 1) - P(m, k), P(m, k + 1) - P(m, k)
        assert math.degrees(math.acos(np.dot(v1, v2) / 2.25)) == pytest.approx(120, abs=1e-3)


def test_2d_fragment_replaces_end_atom_and_stays_flat():
    m = layout.smiles_2d("CC")["molecule"]
    r = edit2d(m, type="fragment", atom=1, key="3")
    t = r["molecule"]
    assert smiles_of(t) == "Cc1ccccc1"
    assert flat(t) and all(abs(x - 1.5) < 0.05 for x in bond_lengths(t))


@pytest.mark.parametrize("frag", FRAGMENTS, ids=lambda f: f.key)
def test_every_fragment_works_in_2d(frag):
    m = layout.smiles_2d("CC")["molecule"]
    r = edit2d(m, type="fragment", atom=1, key=frag.key)
    assert flat(r["molecule"])
    smiles_of(r["molecule"])  # valid structure


def test_2d_ring_fusion_naphthalene():
    m = layout.smiles_2d("c1ccccc1")["molecule"]
    b = next(b for b in m["bonds"] if b["order"] == 1)
    r = edit2d(m, type="fuse", a=b["a"], b=b["b"], size=6, aromatic=True)
    assert smiles_of(r["molecule"]) == "c1ccc2ccccc2c1"
    assert flat(r["molecule"])
    # The new ring is on the outside: no atom closer than a bond length to another non-bonded atom.
    mm = r["molecule"]
    bonded = {frozenset((x["a"], x["b"])) for x in mm["bonds"]}
    for i in range(len(mm["atoms"])):
        for j in range(i + 1, len(mm["atoms"])):
            if frozenset((i, j)) not in bonded:
                assert np.linalg.norm(P(mm, i) - P(mm, j)) > 1.4


def test_clean_2d_keeps_position():
    m = layout.smiles_2d("CCCCO")["molecule"]
    moved = {"atoms": [{**a, "x": a["x"] + 10} for a in m["atoms"]], "bonds": m["bonds"]}
    c = layout.clean_2d(moved)["molecule"]
    cx = sum(a["x"] for a in c["atoms"]) / len(c["atoms"])
    assert cx == pytest.approx(sum(a["x"] for a in moved["atoms"]) / len(moved["atoms"]), abs=1e-6)


def test_2d_export_keeps_implicit_hydrogens_and_cdxml_round_trip():
    m = layout.smiles_2d("CCO")["molecule"]
    block = chem.export_text(m, "mol", dim=2)
    back = Chem.MolFromMolBlock(block)
    assert Chem.MolToSmiles(back) == "CCO"
    assert " 2D" in block.splitlines()[1] or block.splitlines()[1].strip().endswith("2D")
    cdxml = chem.export_text(m, "cdxml", dim=2)
    assert Chem.MolToSmiles(Chem.MolsFromCDXML(cdxml)[0]) == "CCO"
    with pytest.raises(ChemError, match="2D drawing"):
        chem.export_text(m, "xyz", dim=2)


def test_cdxml_import_as_2d_and_3d():
    cdxml = Chem.MolToCDXMLBlock(Chem.MolFromSmiles("c1ccccc1O"))
    two, _ = chem.import_text("phenol.cdxml", cdxml, dim=2)
    assert flat(two) and len(two["atoms"]) == 7
    three, notes = chem.import_text("phenol.cdxml", cdxml, dim=3)
    assert len(three["atoms"]) == 13 and "2D" in notes[0]
    with pytest.raises(ChemError, match="cdxml"):
        chem.import_text("x.cdx", "VjCD0100")
