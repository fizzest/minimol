import stat

import pytest

from minimol import check, chem, gxtb
from minimol.chem import ChemError


def mol(atoms, bonds=()):
    return {
        "atoms": [{"el": el, "x": x, "y": y, "z": z, "charge": 0} for el, x, y, z in atoms],
        "bonds": [{"a": a, "b": b, "order": 1} for a, b in bonds],
    }


def test_a_normal_molecule_has_no_problems():
    water = mol([("O", 0, 0, 0), ("H", 0.96, 0, 0), ("H", -0.24, 0.93, 0)], [(0, 1), (0, 2)])
    assert check.preflight(water, 0) == {"errors": [], "warnings": []}
    ethanol, _ = chem.from_smiles("CCO")
    assert check.preflight(ethanol, 0) == {"errors": [], "warnings": []}


def test_atoms_on_top_of_each_other_are_an_error():
    two_methanes = mol([("C", 0, 0, 0), ("C", 0.2, 0, 0)])
    result = check.preflight(two_methanes, 0)
    assert result["errors"] and "C1–C2" in result["errors"][0] and "0.20" in result["errors"][0]


def test_close_nonbonded_atoms_are_a_clash_with_a_hint():
    # H and C 1.4 A apart without a bond (compare the C-H 1.41 A clash from a bad template).
    clash = mol([("C", 0, 0, 0), ("H", 1.09, 0, 0), ("C", 4, 0, 0), ("H", 5.0, 0, 0), ("H", 2.5, 0.0, 0.0)], [(0, 1), (2, 3)])
    result = check.preflight(clash, 0)
    assert not result["errors"]
    assert any("too close" in w and "UFF or MMFF94" in w for w in result["warnings"])


def test_angle_ends_and_bonded_atoms_are_not_clashes():
    # The two hydrogens of a methylene group are 1.78 A apart: close, but normal (1,3 pair).
    methane, _ = chem.from_smiles("C")
    assert check.preflight(methane, 0)["warnings"] == []


def test_stretched_and_compressed_bonds_are_reported():
    stretched = mol([("C", 0, 0, 0), ("C", 2.4, 0, 0)], [(0, 1)])
    compressed = mol([("C", 0, 0, 0), ("C", 0.9, 0, 0)], [(0, 1)])
    for molecule in (stretched, compressed):
        warnings = check.preflight(molecule, 0)["warnings"]
        assert any("unusual length" in w for w in warnings)


def test_odd_electron_count_is_explained():
    methyl, _ = chem.from_smiles("[CH3]")
    warnings = check.preflight(methyl, 0)["warnings"]
    assert len(warnings) == 1 and "odd (9)" in warnings[0] and "doublet" in warnings[0]
    assert check.preflight(methyl, 0, uhf=1)["warnings"] == []  # open shell on purpose
    assert check.preflight(methyl, 1)["warnings"] == []  # even electron count


def test_odd_electron_count_names_the_atoms_with_open_valences():
    # Benzene with one hydrogen removed: C6H5 radical, the carbon has an open valence.
    phenyl, _ = chem.from_smiles("c1ccccc1")
    h = next(i for i, a in enumerate(phenyl["atoms"]) if a["el"] == "H")
    keep = [i for i in range(len(phenyl["atoms"])) if i != h]
    new = {old: k for k, old in enumerate(keep)}
    radical = {
        "atoms": [phenyl["atoms"][i] for i in keep],
        "bonds": [{**b, "a": new[b["a"]], "b": new[b["b"]]} for b in phenyl["bonds"] if h not in (b["a"], b["b"])],
    }
    warning = check.preflight(radical, 0)["warnings"][0]
    assert "open valences: C" in warning


# ------------------------------------------------- inside the g-xTB call


def test_overlapping_atoms_stop_gxtb_before_it_runs(monkeypatch):
    monkeypatch.setattr(gxtb, "find_binary", lambda refresh=False: pytest.fail("g-xTB must not be started"))
    overlapping = mol([("C", 0, 0, 0), ("C", 0.1, 0, 0), ("H", 5, 0, 0)])
    with pytest.raises(ChemError, match="on top of each other"):
        gxtb.optimize(overlapping)


@pytest.mark.skipif(gxtb.os.name == "nt", reason="uses a shell script as a stand-in for g-xTB")
def test_warnings_come_back_with_the_result_and_with_failures(tmp_path):
    odd = mol([("C", 0, 0, 0), ("H", 1.09, 0, 0), ("H", -0.5, 0.95, 0), ("H", -0.5, -0.95, 0)], [(0, 1), (0, 2), (0, 3)])
    ok = tmp_path / "ok-xtb"
    ok.write_text('#!/bin/sh\nif [ "$1" = "--help" ]; then echo "--gxtb"; exit 0; fi\n'
                  'printf "4\\n energy: -3.5 gnorm: 0.001 xtb: 6.7\\nC 0 0 0\\nH 1.09 0 0\\nH -0.5 0.95 0\\nH -0.5 -0.95 0\\n" > xtbopt.xyz\n'
                  'echo "GEOMETRY OPTIMIZATION CONVERGED"\n')
    bad = tmp_path / "bad-xtb"
    bad.write_text('#!/bin/sh\nif [ "$1" = "--help" ]; then echo "--gxtb"; exit 0; fi\necho "abnormal termination"; exit 1\n')
    for path in (ok, bad):
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
    result = gxtb.optimize(odd, binary=str(ok))
    assert any("electron count is odd" in n for n in result["notes"])
    with pytest.raises(ChemError, match="g-xTB failed.*electron count is odd"):
        gxtb.optimize(odd, binary=str(bad))
