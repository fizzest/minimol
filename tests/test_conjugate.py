import random

import pytest

from minimol import chem, layout
from minimol.chem import ChemError
from minimol.conjugate import max_matching
from minimol.edit import apply_edit

from test_edit import formula


def conjugate(mol, dim=3, **extra):
    return apply_edit(mol, {"type": "conjugate", **extra}, keep_h=True, dim=dim)


def doubles(mol):
    return sum(1 for b in mol["bonds"] if b["order"] == 2)


def drawing(smiles):
    return layout.smiles_2d(smiles)["molecule"]


def single_bonds(mol):
    return {**mol, "bonds": [{**b, "order": 1} for b in mol["bonds"]]}


def hydrogens_2d(mol):
    return sum(a.get("h", 0) for a in mol["atoms"])


# ------------------------------------------------------------ the matching


def brute_force_size(edges):
    best = 0

    def rec(start, used, count):
        nonlocal best
        best = max(best, count)
        for k in range(start, len(edges)):
            a, b = edges[k]
            if a not in used and b not in used:
                rec(k + 1, used | {a, b}, count + 1)

    rec(0, frozenset(), 0)
    return best


def test_matching_is_maximum_on_random_graphs_with_odd_cycles():
    rng = random.Random(7)
    for _ in range(400):
        n = rng.randint(2, 10)
        edges = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < rng.choice([0.2, 0.35, 0.5])]
        start = []
        used = set()
        for a, b in rng.sample(edges, len(edges)):
            if a not in used and b not in used and rng.random() < 0.5:
                start.append((a, b))
                used |= {a, b}
        mate = max_matching(n, edges, start)
        pairs = {(min(v, mate[v]), max(v, mate[v])) for v in range(n) if mate[v] != -1}
        assert all(mate[mate[v]] == v for v in range(n) if mate[v] != -1)
        assert pairs <= {(min(a, b), max(a, b)) for a, b in edges}
        assert len(pairs) == brute_force_size(edges)
        assert all(mate[a] != -1 and mate[b] != -1 for a, b in start)  # earlier pairs are never given up


# ------------------------------------------------------------ the edit


def test_cyclohexane_becomes_benzene_in_3d_and_2d():
    cyclohexane, _ = chem.from_smiles("C1CCCCC1")
    r = conjugate(cyclohexane)
    assert formula(r["molecule"]) == "C6H6" and doubles(r["molecule"]) == 3
    assert any("Removed 6 hydrogen" in n for n in r["notes"])
    r2 = conjugate(drawing("C1CCCCC1"), dim=2)
    assert doubles(r2["molecule"]) == 3 and hydrogens_2d(r2["molecule"]) == 6


def test_without_ring_hydrogen_removal_nothing_dehydrogenates():
    cyclohexane, _ = chem.from_smiles("C1CCCCC1")
    r = conjugate(cyclohexane, dehydro=False)
    assert formula(r["molecule"]) == "C6H12" and doubles(r["molecule"]) == 0
    assert r["notes"] == ["Already as conjugated as it can be."]


def test_tetralin_becomes_naphthalene():
    tetralin, _ = chem.from_smiles("c1ccc2CCCCc2c1")
    r = conjugate(tetralin)
    assert formula(r["molecule"]) == "C10H8" and doubles(r["molecule"]) == 5


def test_things_that_must_not_change_stay_the_same():
    for smiles, expected in [
        ("CCc1ccccc1", "C8H10"),  # ethyl group: not in a ring
        ("O=C(C)c1ccccc1", "C8H8O"),  # carbonyl
        ("c1cc[nH]c1", "C4H5N"),  # N-H of pyrrole
        ("C#Cc1ccccc1", "C8H6"),  # triple bond
        ("Oc1ccccc1", "C6H6O"),  # O-H
        ("CC=CC=CC", "C6H10"),  # open chain: nothing to aromatize
    ]:
        m, _ = chem.from_smiles(smiles)
        r = conjugate(m)
        assert formula(r["molecule"]) == expected, smiles
        assert r["notes"] == ["Already as conjugated as it can be."], smiles


def test_a_broken_kekule_structure_is_repaired_without_changing_hydrogens():
    naphthalene, _ = chem.from_smiles("c1ccc2ccccc2c1")
    first_double = next(i for i, b in enumerate(naphthalene["bonds"]) if b["order"] == 2)
    broken = {**naphthalene, "bonds": [{**b, "order": 2 if i == first_double else 1} for i, b in enumerate(naphthalene["bonds"])]}
    r = conjugate(broken)
    assert formula(r["molecule"]) == "C10H8" and doubles(r["molecule"]) == 5


def test_large_pah_with_every_bond_single_gets_a_perfect_alternation():
    coronene = drawing("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67")
    r = conjugate(single_bonds(coronene), dim=2)
    assert doubles(r["molecule"]) == 12 and hydrogens_2d(r["molecule"]) == 12


def test_only_the_chosen_atoms_change():
    two_rings, _ = chem.from_smiles("C1CCCCC1C1CCCCC1")
    first_ring = list(range(6))
    r = conjugate(two_rings, atoms=first_ring)
    assert formula(r["molecule"]) == "C12H16"  # one benzene ring, one cyclohexyl
    assert doubles(r["molecule"]) == 3


def test_ring_carbon_without_partner_stays_and_is_not_called_a_radical():
    acridane = drawing("C1=CC=C2C(=C1)CC1=CC=CC=C1N2")
    r = conjugate(acridane, dim=2)
    assert doubles(r["molecule"]) == 6 and hydrogens_2d(r["molecule"]) == 11  # C13H11N: CH2 and N-H stay
    assert r["notes"] == ["Already as conjugated as it can be."]


def test_open_valence_centres_are_reported():
    # Allyl radical: three carbons cannot all take a double bond, so one stays open-valence.
    allyl, _ = chem.from_smiles("[CH2]C=C")
    r = conjugate(allyl)
    assert formula(r["molecule"]) == "C3H5" and doubles(r["molecule"]) == 1
    assert any("radical" in n and "C1" in n for n in r["notes"])


def test_bad_arguments():
    water, _ = chem.from_smiles("O")
    with pytest.raises(ChemError, match="dehydro"):
        conjugate(water, dehydro="yes")
    with pytest.raises(ChemError, match="list"):
        conjugate(water, atoms=3)
    with pytest.raises(ChemError, match="does not exist"):
        conjugate(water, atoms=[99])
    with pytest.raises(ChemError, match="no non-hydrogen"):
        conjugate({"atoms": [{"el": "H", "x": 0, "y": 0, "z": 0, "charge": 0}], "bonds": []})


def test_result_is_a_valid_molecule_with_unchanged_geometry_of_untouched_atoms():
    m, _ = chem.from_smiles("CCc1ccccc1C1CCCCC1")
    r = conjugate(m)
    chem.validate_json(r["molecule"])
    for old, new in enumerate(r["map"]):
        if new >= 0 and old < 2:  # the ethyl carbons were not touched
            assert r["molecule"]["atoms"][new]["x"] == m["atoms"][old]["x"]
