import math
import os
import stat
import textwrap

import pytest

from minimol import chem, gxtb
from minimol.chem import ChemError

BINARY = gxtb.find_binary(refresh=True)
needs_gxtb = pytest.mark.skipif(BINARY is None, reason="g-xTB binary not installed")
# Some tests stand in for g-xTB with small shell scripts, which Windows cannot run.
needs_sh = pytest.mark.skipif(os.name == "nt", reason="uses a shell script as a stand-in for g-xTB")

WATER = {
    "atoms": [
        {"el": "O", "x": 0.0, "y": 0.0, "z": 0.0, "charge": 0},
        {"el": "H", "x": 1.10, "y": 0.0, "z": 0.0, "charge": 0},
        {"el": "H", "x": -0.2, "y": 1.05, "z": 0.1, "charge": 0},
    ],
    "bonds": [{"a": 0, "b": 1, "order": 1}, {"a": 0, "b": 2, "order": 1}],
}


def d(mol, i, j):
    a, b = mol["atoms"][i], mol["atoms"][j]
    return math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))


# ------------------------------------------------------ binary detection


def _fake_binary(tmp_path, name, help_text, script=""):
    path = tmp_path / name
    path.write_text(
        "#!/bin/sh\n"
        f'if [ "$1" = "--help" ]; then printf %s "{help_text}"; exit 0; fi\n' + script
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


@needs_sh
def test_env_var_binary_must_support_gxtb(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    plain = _fake_binary(tmp_path, "plain-xtb", "--gfn INT")
    monkeypatch.setenv(gxtb.ENV_VAR, plain)
    assert gxtb.find_binary(refresh=True) is None
    good = _fake_binary(tmp_path, "my-gxtb", "--gxtb  use the g-xTB method")
    monkeypatch.setenv(gxtb.ENV_VAR, good)
    assert gxtb.find_binary(refresh=True) == good
    gxtb.find_binary(refresh=True)  # leave the cache as the fixture found it


@pytest.fixture(autouse=True)
def _restore_cache():
    yield
    gxtb._cache.clear()


def test_missing_binary_message(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.delenv(gxtb.ENV_VAR, raising=False)
    monkeypatch.setattr(gxtb.gxtb_install, "installed_binary", lambda: None)  # nor a copy installed by the app
    gxtb._cache.clear()
    with pytest.raises(ChemError, match="g-xTB was not found"):
        gxtb.optimize(WATER)


# ------------------------------------------------------- output parsing


def test_error_message_puts_innermost_cause_first():
    log = textwrap.dedent(
        """\
        -2- reading geometry input 'in.xyz' failed
        -1- Error: Atom number missmatch in xyz file
        """
    )
    assert gxtb._error_message(log, 128) == (
        "Error: Atom number missmatch in xyz file; reading geometry input 'in.xyz' failed"
    )
    assert "exit status 3" in gxtb._error_message("something\nbad happened\n", 3)


def test_read_xtbopt_checks_atoms():
    text = "3\n energy: -76.4 gnorm: 0.0003 xtb: 6.7.1\nO 0 0 0\nH 0.96 0 0\nH 0 0.96 0\n"
    coords, energy, gnorm = gxtb._read_xtbopt(text, WATER["atoms"])
    assert energy == -76.4 and gnorm == 0.0003 and coords[1] == (0.96, 0.0, 0.0)
    with pytest.raises(ChemError, match="element of atom 2"):
        gxtb._read_xtbopt(text.replace("H 0.96", "F 0.96"), WATER["atoms"])
    with pytest.raises(ChemError, match="number of atoms"):
        gxtb._read_xtbopt("2\n energy: 1\nO 0 0 0\nH 1 0 0\n", WATER["atoms"])


@needs_sh
def test_fake_binary_full_run(tmp_path, monkeypatch):
    """Exercise the file handling with a stand-in that mimics xtb's outputs."""
    script = textwrap.dedent(
        """\
        printf '3\\n energy: -1.5 gnorm: 0.0001 xtb: fake\\nO 0 0 0\\nH 0.95 0 0\\nH 0 0.95 0\\n' > xtbopt.xyz
        echo ':: total energy  -1.25 Eh ::'
        echo '*** GEOMETRY OPTIMIZATION CONVERGED AFTER 3 ITERATIONS ***'
        echo "ARGS $@"
        """
    )
    fake = _fake_binary(tmp_path, "gxtb", "--gxtb", script)
    result = gxtb.optimize(WATER, binary=fake, uhf=0, fixed=[0])
    assert result["converged"] and result["energy"] == -1.5 and result["initial_energy"] == -1.25
    assert result["energy_unit"] == "Eh" and result["method"] == "g-xTB"
    assert result["molecule"]["bonds"] == WATER["bonds"]
    assert result["molecule"]["atoms"][1]["x"] == 0.95


@needs_sh
def test_fake_binary_failure_is_reported(tmp_path):
    fake = _fake_binary(tmp_path, "gxtb", "--gxtb", "echo '-1- Error: SCF did not converge'\nexit 1\n")
    with pytest.raises(ChemError, match="g-xTB failed: Error: SCF did not converge"):
        gxtb.optimize(WATER, binary=fake)


# ----------------------------------------------------- input validation


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"uhf": 1}, "impossible"),
        ({"charge": -1, "uhf": 0}, "impossible"),
        ({"charge": 1.5}, "integer"),
        ({"level": "fast"}, "level"),
        ({"fixed": [0, 1, 2]}, "All atoms are fixed"),
        ({"fixed": [7]}, "out of range"),
        ({"max_cycles": 0}, "max_cycles"),
    ],
)
def test_validation(kwargs, message):
    with pytest.raises(ChemError, match=message):
        gxtb.optimize(WATER, binary="/bin/false", **kwargs)


def test_elements_beyond_lr_are_rejected():
    rf = {"atoms": [{"el": "Rf", "x": 0, "y": 0, "z": 0, "charge": 0}], "bonds": []}
    with pytest.raises(ChemError, match="atom 1 is Rf"):
        gxtb.optimize(rf, binary="/bin/false")


def test_periodic_table_matches_rdkit():
    from rdkit import Chem

    pt = Chem.GetPeriodicTable()
    assert len(gxtb._PERIODIC) == 103
    assert all(pt.GetElementSymbol(i + 1) == s for i, s in enumerate(gxtb._PERIODIC))


# ------------------------------------------------- real g-xTB (optional)


@needs_gxtb
def test_real_water_optimization():
    result = gxtb.optimize(WATER)
    mol = result["molecule"]
    assert result["converged"]
    assert result["energy"] < result["initial_energy"] < 0
    assert 0.93 < d(mol, 0, 1) < 0.99 and 0.93 < d(mol, 0, 2) < 0.99


@needs_gxtb
def test_real_fixed_atoms_and_charge():
    result = gxtb.optimize(WATER, fixed=[0, 1])
    assert result["molecule"]["atoms"][1]["x"] == pytest.approx(1.10, abs=1e-6)
    hydroxide = {"atoms": WATER["atoms"][:2], "bonds": WATER["bonds"][:1]}
    hydroxide["atoms"] = [dict(hydroxide["atoms"][0], charge=-1), hydroxide["atoms"][1]]
    r = gxtb.optimize(hydroxide)
    assert r["charge"] == -1 and r["converged"]


@needs_gxtb
def test_real_non_convergence_is_reported():
    result = gxtb.optimize(WATER, max_cycles=1)
    assert not result["converged"] and "did not converge" in result["notes"][0]


@needs_gxtb
def test_real_smiles_structure_optimizes():
    mol, _ = chem.from_smiles("CC(=O)O")
    result = gxtb.optimize(mol)
    assert result["converged"]
    co = [b for b in result["molecule"]["bonds"] if b["order"] == 2][0]
    assert 1.18 < d(result["molecule"], co["a"], co["b"]) < 1.24  # C=O
