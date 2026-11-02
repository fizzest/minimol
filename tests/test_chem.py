import math

import pytest
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers as FFH

from minimol import chem
from minimol.chem import ChemError


def atom(el, x=0.0, y=0.0, z=0.0, charge=0):
    return {"el": el, "x": x, "y": y, "z": z, "charge": charge}


def dist(mol, i, j):
    a, b = mol["atoms"][i], mol["atoms"][j]
    return math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))


ETHANE_ROUGH = {
    "atoms": [atom("C"), atom("C", 1.9, 0.2)],
    "bonds": [{"a": 0, "b": 1, "order": 1}],
}


# ------------------------------------------------------------ validation


@pytest.mark.parametrize(
    "bad, message",
    [
        ({"atoms": [atom("Xx")], "bonds": []}, "Unknown element"),
        ({"atoms": [atom("c")], "bonds": []}, "Invalid element"),
        ({"atoms": [{**atom("C"), "x": float("nan")}], "bonds": []}, "finite"),
        ({"atoms": [{**atom("C"), "x": "1"}], "bonds": []}, "number"),
        ({"atoms": [atom("C")], "bonds": [{"a": 0, "b": 1, "order": 1}]}, "missing atom"),
        ({"atoms": [atom("C")], "bonds": [{"a": 0, "b": 0, "order": 1}]}, "itself"),
        ({"atoms": [atom("C"), atom("C", 1.5)], "bonds": [{"a": 0, "b": 1, "order": 4}]}, "order"),
        (
            {"atoms": [atom("C"), atom("C", 1.5)], "bonds": [{"a": 0, "b": 1, "order": 1}, {"a": 1, "b": 0, "order": 1}]},
            "twice",
        ),
        ({"atoms": [{**atom("C"), "charge": 1.5}], "bonds": []}, "integer"),
        ([], "object"),
    ],
)
def test_validate_rejects_bad_input(bad, message):
    with pytest.raises(ChemError, match=message):
        chem.validate_json(bad)


def test_round_trip_keeps_atoms_bonds_and_charges():
    data = {
        "atoms": [atom("N", charge=1), atom("H", 1.0), atom("H", -0.3, 0.95), atom("H", -0.3, -0.5, 0.8), atom("H", -0.3, -0.5, -0.8)],
        "bonds": [{"a": 0, "b": i, "order": 1} for i in range(1, 5)],
    }
    out = chem.mol_to_json(chem.mol_from_json(data))
    assert out == data


def test_wysiwyg_open_valence_is_a_radical_and_is_reported():
    data = {"atoms": [atom("C"), atom("H", 1.09), atom("H", -0.5, 0.9), atom("H", -0.5, -0.9)], "bonds": [{"a": 0, "b": i, "order": 1} for i in (1, 2, 3)]}
    mol = chem.mol_from_json(data)
    assert mol.GetAtomWithIdx(0).GetNumRadicalElectrons() == 1
    result = chem.optimize(data)
    assert any("open valences" in n and "atoms 1" in n for n in result["notes"])


def test_pentavalent_carbon_reports_one_based_atom_number():
    data = {"atoms": [atom("H"), atom("C", 1.0)] + [atom("H", 1.0, 1.0 * k) for k in range(1, 5)], "bonds": [{"a": 1, "b": i, "order": 1} for i in (0, 2, 3, 4, 5)]}
    with pytest.raises(ChemError, match=r"atom 2 C"):
        chem.mol_from_json(data)


# ------------------------------------------------------------ hydrogens


def test_add_hydrogens_fills_valences_and_keeps_indices():
    data = {"atoms": [atom("C"), atom("O", 1.43)], "bonds": [{"a": 0, "b": 1, "order": 1}]}
    out, notes = chem.add_hydrogens(data)
    els = [a["el"] for a in out["atoms"]]
    assert els[:2] == ["C", "O"]
    assert sorted(els[2:]) == ["H"] * 4
    assert out["atoms"][0]["x"] == 0 and out["atoms"][1]["x"] == 1.43
    assert "Added 4" in notes[0]
    # Hydrogens sit at bonding distance from their heavy atom.
    for b in out["bonds"]:
        if out["atoms"][b["b"]]["el"] == "H":
            assert 0.9 < dist(out, b["a"], b["b"]) < 1.2


def test_add_hydrogens_respects_formal_charge_and_existing_h():
    data = {"atoms": [atom("N", charge=1), atom("H", 1.0)], "bonds": [{"a": 0, "b": 1, "order": 1}]}
    out, _ = chem.add_hydrogens(data)
    assert [a["el"] for a in out["atoms"]].count("H") == 4  # NH4+


def test_add_hydrogens_on_saturated_molecule_adds_nothing():
    methane, _ = chem.from_smiles("C")
    out, notes = chem.add_hydrogens(methane)
    assert len(out["atoms"]) == 5 and "No open valences" in notes[0]


# ----------------------------------------------------------- optimization


@pytest.mark.parametrize("method", ["MMFF94", "MMFF94s", "UFF"])
def test_optimize_ethane_gives_sane_geometry(method):
    rough, _ = chem.add_hydrogens(ETHANE_ROUGH)
    result = chem.optimize(rough, method)
    mol = result["molecule"]
    assert result["converged"]
    assert result["energy"] < result["initial_energy"]
    assert 1.49 < dist(mol, 0, 1) < 1.56
    assert [a["el"] for a in mol["atoms"]] == [a["el"] for a in rough["atoms"]]
    assert mol["bonds"] == rough["bonds"]


def test_mmff_energy_matches_rdkit_directly():
    smi, _ = chem.from_smiles("CC(=O)Nc1ccc(O)cc1")
    result = chem.optimize(smi, "MMFF94")
    mol = chem.mol_from_json(result["molecule"])
    props = FFH.MMFFGetMoleculeProperties(mol)
    ff = FFH.MMFFGetMoleculeForceField(mol, props, ignoreInterfragInteractions=False)
    assert ff.CalcEnergy() == pytest.approx(result["energy"], abs=1e-3)


def test_fragments_interact_during_optimization():
    # Two waters 1.2 A apart (O...H) would collapse into each other if
    # inter-fragment interactions were ignored, or be pushed apart if not.
    water, _ = chem.from_smiles("O")
    two = {"atoms": water["atoms"] + [{**a, "x": a["x"] + 1.6} for a in water["atoms"]], "bonds": water["bonds"] + [{"a": b["a"] + 3, "b": b["b"] + 3, "order": 1} for b in water["bonds"]]}
    result = chem.optimize(two, "MMFF94")
    o1, o2 = result["molecule"]["atoms"][0], result["molecule"]["atoms"][3]
    assert math.dist((o1["x"], o1["y"], o1["z"]), (o2["x"], o2["y"], o2["z"])) > 2.4


def test_fixed_atoms_do_not_move():
    rough, _ = chem.add_hydrogens(ETHANE_ROUGH)
    result = chem.optimize(rough, "MMFF94", fixed=[0, 1])
    for i in (0, 1):
        for k in "xyz":
            assert result["molecule"]["atoms"][i][k] == pytest.approx(rough["atoms"][i][k], abs=1e-6)


def test_optimize_errors():
    with pytest.raises(ChemError, match="empty"):
        chem.optimize({"atoms": [], "bonds": []})
    with pytest.raises(ChemError, match="Unknown force field"):
        chem.optimize(ETHANE_ROUGH, "AMBER")
    with pytest.raises(ChemError, match="out of range"):
        chem.optimize(ETHANE_ROUGH, fixed=[5])
    uranium = {"atoms": [atom("U"), atom("F", 2.0)], "bonds": [{"a": 0, "b": 1, "order": 1}]}
    with pytest.raises(ChemError, match="Try UFF"):
        chem.optimize(uranium, "MMFF94")


# ------------------------------------------------------------------ SMILES


def test_from_smiles_builds_3d_with_hydrogens():
    mol, _ = chem.from_smiles("c1ccccc1")
    assert len(mol["atoms"]) == 12
    orders = sorted(b["order"] for b in mol["bonds"])
    assert orders == [1] * 9 + [2] * 3  # Kekulé form
    assert any(abs(a["z"]) > 1e-3 for a in mol["atoms"]) or any(abs(a["y"]) > 1e-3 for a in mol["atoms"])


@pytest.mark.parametrize("bad", ["", "   ", "C1CC", "not a smiles(", 5])
def test_from_smiles_rejects_bad_input(bad):
    with pytest.raises(ChemError):
        chem.from_smiles(bad)


# ---------------------------------------------------------- import/export

BENZENE_XYZ = """12
benzene
C 0.000 1.396 0.000
C 1.209 0.698 0.000
C 1.209 -0.698 0.000
C 0.000 -1.396 0.000
C -1.209 -0.698 0.000
C -1.209 0.698 0.000
H 0.000 2.480 0.000
H 2.148 1.240 0.000
H 2.148 -1.240 0.000
H 0.000 -2.480 0.000
H -2.148 -1.240 0.000
H -2.148 1.240 0.000
"""


def test_import_xyz_perceives_bond_orders():
    mol, notes = chem.import_text("benzene.xyz", BENZENE_XYZ)
    assert notes == []
    assert len(mol["atoms"]) == 12 and len(mol["bonds"]) == 12
    assert sorted(b["order"] for b in mol["bonds"]).count(2) == 3


def test_import_xyz_tolerates_case_atomic_numbers_and_extra_columns():
    text = "3\nwater\n8 0.0 0.0 0.0 extra\nh 0.96 0.0 0.0\nH -0.24 0.93 0.0 0.1\n"
    mol, _ = chem.import_text("w.XYZ", text)
    assert [a["el"] for a in mol["atoms"]] == ["O", "H", "H"]
    assert len(mol["bonds"]) == 2


def test_import_xyz_multiframe_loads_first_and_says_so():
    frame = "1\nf\nHe 0 0 0\n"
    mol, notes = chem.import_text("t.xyz", frame + frame.replace("0 0 0", "1 1 1"))
    assert mol["atoms"][0]["x"] == 0 and "first" in notes[0]


def test_import_xyz_falls_back_to_connectivity_when_bond_orders_fail(monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("no valid bond order assignment")

    monkeypatch.setattr(chem.rdDetermineBonds, "DetermineBonds", fail)
    mol, notes = chem.import_text("benzene.xyz", BENZENE_XYZ)
    assert len(mol["bonds"]) == 12
    assert all(b["order"] == 1 for b in mol["bonds"])
    assert any("set to single" in n for n in notes)


@pytest.mark.parametrize(
    "name, text, message",
    [
        ("a.xyz", "x\n", "number of atoms"),
        ("a.xyz", "3\n\nC 0 0 0\n", "shorter"),
        ("a.xyz", "1\n\nC 0 0 zz\n", "not a number"),
        ("a.xyz", "1\n\nQq 0 0 0\n", "Unknown element"),
        ("a.txt", "1\n\nC 0 0 0\n", "Unsupported file type"),
        ("noext", "x", "extension"),
        ("a.xyz", "   ", "empty"),
        ("a.mol", "garbage\n\n\n", "MOL/SDF"),
    ],
)
def test_import_errors(name, text, message):
    with pytest.raises(ChemError, match=message):
        chem.import_text(name, text)


def test_sdf_roundtrip_and_multi_record_note():
    mol, _ = chem.from_smiles("CC(=O)O")
    sdf = chem.export_text(mol, "sdf", "acetic acid")
    assert sdf.startswith("acetic acid\n") and sdf.rstrip().endswith("$$$$")
    back, notes = chem.import_text("x.sdf", sdf + sdf)
    assert back["bonds"] == mol["bonds"]
    for a, b in zip(back["atoms"], mol["atoms"]):
        assert a["el"] == b["el"] and a["charge"] == b["charge"]
        # MOL/SDF (V2000) stores coordinates with 4 decimals.
        assert [a[k] for k in "xyz"] == pytest.approx([b[k] for k in "xyz"], abs=1e-4)
    assert "2 molecules" in notes[0]


def test_2d_molfile_is_embedded_in_3d():
    block = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))  # 2D, no H, all z = 0
    mol, notes = chem.import_text("e.mol", block)
    assert len(mol["atoms"]) == 9
    assert any(abs(a["z"]) > 0.1 for a in mol["atoms"])
    assert "2D" in notes[0]


def test_xyz_export_import_roundtrip():
    mol, _ = chem.from_smiles("C#N")
    text = chem.export_text(mol, "xyz")
    back, _ = chem.import_text("h.xyz", text)
    assert [a["el"] for a in back["atoms"]] == [a["el"] for a in mol["atoms"]]
    assert sorted(b["order"] for b in back["bonds"]) == [1, 3]


def test_pdb_and_mol2_import():
    mol, _ = chem.from_smiles("CCO")
    rd = chem.mol_from_json(mol)
    pdb_mol, _ = chem.import_text("e.pdb", Chem.MolToPDBBlock(rd))
    assert len(pdb_mol["atoms"]) == 9
    mol2 = """@<TRIPOS>MOLECULE
water
 3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 O1          0.0000    0.0000    0.0000 O.3     1  HOH1        0.0000
      2 H1          0.9572    0.0000    0.0000 H       1  HOH1        0.0000
      3 H2         -0.2400    0.9266    0.0000 H       1  HOH1        0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
"""
    water, _ = chem.import_text("w.mol2", mol2)
    assert [a["el"] for a in water["atoms"]] == ["O", "H", "H"]


def test_export_of_chemically_invalid_structure_still_works():
    bad = {"atoms": [atom("C")] + [atom("H", 1.0, k) for k in range(5)], "bonds": [{"a": 0, "b": i, "order": 1} for i in range(1, 6)]}
    assert "C" in chem.export_text(bad, "mol")
    assert chem.export_text(bad, "xyz").startswith("6\n")


def test_export_rejects_unknown_format():
    with pytest.raises(ChemError, match="format"):
        chem.export_text(ETHANE_ROUGH, "cif")
