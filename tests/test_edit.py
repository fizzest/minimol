import math
from collections import Counter

import numpy as np
import pytest

from minimol import chem
from minimol.chem import ChemError
from minimol.edit import FRAGMENTS, apply_edit, template


def formula(mol):
    c = Counter(a["el"] for a in mol["atoms"])
    order = ["C", "H"] + sorted(k for k in c if k not in ("C", "H"))
    return "".join(k + (str(c[k]) if c[k] > 1 else "") for k in order if c[k])


def P(mol, i):
    a = mol["atoms"][i]
    return np.array([a["x"], a["y"], a["z"]])


def d(mol, i, j):
    return float(np.linalg.norm(P(mol, i) - P(mol, j)))


def bonded(mol, i, j):
    return any({b["a"], b["b"]} == {i, j} for b in mol["bonds"])


def neighbors(mol, i):
    return [b["b"] if b["a"] == i else b["a"] for b in mol["bonds"] if i in (b["a"], b["b"])]


def edit(mol, keep_h=True, **op):
    return apply_edit(mol, op, keep_h=keep_h)


def min_nonbonded(mol):
    best = math.inf
    n = len(mol["atoms"])
    for i in range(n):
        for j in range(i + 1, n):
            if not bonded(mol, i, j):
                best = min(best, d(mol, i, j))
    return best


def sane(mol):
    """Valid JSON, passes RDKit sanitization, no atoms on top of each other."""
    chem.validate_json(mol)
    chem.mol_from_json(mol)  # raises on valence problems
    assert min_nonbonded(mol) > 0.9


EMPTY = {"atoms": [], "bonds": []}


def methane():
    return edit(EMPTY, type="add_atom", element="C", pos=[0, 0, 0])["molecule"]


def ethane():
    """Ethane and the index of its second carbon (atom 1 is a hydrogen)."""
    r = edit(methane(), type="grow", atom=0, element="C")
    return r["molecule"], r["focus"]


# ------------------------------------------------------------- new atoms


def test_add_atom_with_hydrogens_is_tetrahedral():
    r = edit(EMPTY, type="add_atom", element="C", pos=[1, 2, 3])
    m = r["molecule"]
    assert formula(m) == "CH4" and r["focus"] == 0
    for h in range(1, 5):
        assert d(m, 0, h) == pytest.approx(0.76 + 0.31, abs=1e-6)
    for i in range(1, 5):
        for j in range(i + 1, 5):
            v1, v2 = P(m, i) - P(m, 0), P(m, j) - P(m, 0)
            angle = math.degrees(math.acos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
            assert angle == pytest.approx(109.47, abs=2)


@pytest.mark.parametrize("el, charge_op, expected", [("N", None, "H3N"), ("O", None, "H2O"), ("N", 1, "H4N"), ("H", None, "H")])
def test_add_atom_other_elements(el, charge_op, expected):
    m = edit(EMPTY, type="add_atom", element=el, pos=[0, 0, 0])["molecule"]
    if charge_op:
        m = edit(m, type="charge", atoms=[0], delta=charge_op)["molecule"]
    assert formula(m) == expected
    sane(m)


def test_add_atom_without_hydrogens():
    m = edit(EMPTY, keep_h=False, type="add_atom", element="C", pos=[0, 0, 0])["molecule"]
    assert formula(m) == "C"


# ------------------------------------------------------- growing / bonding


def test_grow_replaces_a_hydrogen_and_hydrogenates_the_new_atom():
    m = methane()
    r = edit(m, type="grow", atom=0, element="C")
    e = r["molecule"]
    assert formula(e) == "C2H6"
    assert r["map"][0] == 0 and r["map"].count(-1) == 1  # exactly one H was replaced
    j = r["focus"]
    assert e["atoms"][j]["el"] == "C" and bonded(e, 0, j)
    assert d(e, 0, j) == pytest.approx(1.52, abs=1e-6)
    sane(e)


def test_grow_into_an_open_valence_keeps_existing_hydrogens():
    m = {"atoms": [{"el": "C", "x": 0, "y": 0, "z": 0, "charge": 0}], "bonds": []}
    r = edit(m, keep_h=False, type="grow", atom=0, element="O")
    assert formula(r["molecule"]) == "CO" and r["map"] == [0]


def test_grow_follows_the_drag_direction_hint():
    m = methane()
    r = edit(m, type="grow", atom=0, element="C", hint=[0, 0, 1])
    e = r["molecule"]
    u = (P(e, r["focus"]) - P(e, 0)) / 1.52
    assert np.allclose(u, [0, 0, 1])  # placed exactly along the drag direction
    # ...and the hydrogen that pointed most along +z is the one that was replaced.
    removed = r["map"].index(-1)
    assert m["atoms"][removed]["z"] == max(a["z"] for a in m["atoms"][1:])


def test_chain_grows_from_focus():
    m = methane()
    focus = 0
    for _ in range(3):
        r = edit(m, type="grow", atom=focus, element="C")
        m, focus = r["molecule"], r["focus"]
    assert formula(m) == "C4H10"
    sane(m)


def test_grow_on_bonded_hydrogen_is_refused():
    m = methane()
    with pytest.raises(ChemError, match="hydrogen"):
        edit(m, type="grow", atom=1, element="C")


def test_bond_between_two_methanes_removes_facing_hydrogens():
    a = methane()
    b = edit(a, type="add_atom", element="C", pos=[3.0, 0, 0])["molecule"]
    r = edit(b, type="bond", a=0, b=5)
    m = r["molecule"]
    assert formula(m) == "C2H6"
    new_c = r["map"][5]
    assert bonded(m, 0, new_c)
    # The hydrogens that were removed pointed at the other carbon.
    for h in neighbors(m, 0):
        if m["atoms"][h]["el"] == "H":
            assert np.dot(P(m, h) - P(m, 0), [1, 0, 0]) < 0.9


def test_bond_errors():
    m, c2 = ethane()
    with pytest.raises(ChemError, match="already bonded"):
        edit(m, type="bond", a=0, b=c2)
    with pytest.raises(ChemError, match="itself"):
        edit(m, type="bond", a=0, b=0)


# ------------------------------------------------ element, order, charge


def test_element_change_fixes_hydrogen_count():
    m = methane()
    assert formula(edit(m, type="element", atom=0, element="O")["molecule"]) == "H2O"
    assert formula(edit(m, type="element", atom=0, element="N")["molecule"]) == "H3N"
    # Clicking a hydrogen with carbon turns it into a methyl group.
    e = edit(m, type="element", atom=1, element="C")["molecule"]
    assert formula(e) == "C2H6"
    sane(e)


def test_element_change_without_keep_h_only_changes_the_element():
    m = methane()
    assert formula(edit(m, keep_h=False, type="element", atom=0, element="O")["molecule"]) == "H4O"


def test_bond_order_cycle_adjusts_hydrogens():
    m, c2 = ethane()
    seen = []
    for order in (2, 3, 1):
        r = edit(m, type="bond_order", a=0, b=c2, order=order)
        m, c2 = r["molecule"], r["map"][c2]
        seen.append(formula(m))
        sane(m)
    assert seen == ["C2H4", "C2H2", "C2H6"]
    h = next(j for j in neighbors(m, 0) if m["atoms"][j]["el"] == "H")
    with pytest.raises(ChemError, match="single bond"):
        edit(m, type="bond_order", a=0, b=h, order=2)


def test_charge_adjusts_hydrogens():
    nh3 = edit(EMPTY, type="add_atom", element="N", pos=[0, 0, 0])["molecule"]
    nh4 = edit(nh3, type="charge", atoms=[0], delta=1)["molecule"]
    assert formula(nh4) == "H4N" and nh4["atoms"][0]["charge"] == 1
    water = edit(EMPTY, type="add_atom", element="O", pos=[0, 0, 0])["molecule"]
    oh = edit(water, type="charge", atoms=[0], delta=-1)["molecule"]
    assert formula(oh) == "HO"


# -------------------------------------------------------------- deleting


def test_delete_heavy_atom_caps_neighbour_with_hydrogen():
    ethanol, _ = chem.from_smiles("CCO")
    o = next(i for i, a in enumerate(ethanol["atoms"]) if a["el"] == "O")
    r = edit(ethanol, type="delete", atoms=[o])
    assert formula(r["molecule"]) == "C2H6"
    assert r["map"][o] == -1
    sane(r["molecule"])


def test_delete_hydrogen_leaves_an_open_valence():
    m = methane()
    r = edit(m, type="delete", atoms=[1])
    assert formula(r["molecule"]) == "CH3"


def test_delete_bond_caps_both_ends():
    m, c2 = ethane()
    r = edit(m, type="delete", bonds=[[0, c2]])
    assert formula(r["molecule"]) == "C2H8"
    assert len(r["molecule"]["bonds"]) == 8


def test_delete_without_keep_h():
    ethanol, _ = chem.from_smiles("CCO")
    o = next(i for i, a in enumerate(ethanol["atoms"]) if a["el"] == "O")
    assert formula(edit(ethanol, keep_h=False, type="delete", atoms=[o])["molecule"]) == "C2H5"


def test_map_tracks_surviving_atoms():
    ethanol, _ = chem.from_smiles("CCO")
    o = next(i for i, a in enumerate(ethanol["atoms"]) if a["el"] == "O")
    r = edit(ethanol, type="delete", atoms=[o])
    for old, new in enumerate(r["map"]):
        if new >= 0:
            assert r["molecule"]["atoms"][new]["el"] == ethanol["atoms"][old]["el"]
            if ethanol["atoms"][old]["el"] != "H":
                assert np.allclose(P(r["molecule"], new), P(ethanol, old))


# -------------------------------------------------------------- fragments


def test_3_on_a_methyl_turns_it_into_phenyl():
    # The user's case: hover the CH3 of ethane, press 3 -> toluene (not benzyl-methyl).
    m, c2 = ethane()
    r = edit(m, type="fragment", atom=c2, key="3")
    t = r["molecule"]
    assert formula(t) == "C7H8"
    sane(t)
    # The target moves to the para carbon so key sequences continue from there.
    para = r["focus"]
    assert t["atoms"][para]["el"] == "C"
    assert sum(1 for j in neighbors(t, para) if t["atoms"][j]["el"] == "H") == 1
    assert d(t, 0, para) > 4.0  # far side of the ring from the original carbon


def test_3_on_methane_gives_benzene():
    assert formula(edit(methane(), type="fragment", atom=0, key="3")["molecule"]) == "C6H6"


def test_1_adds_a_carbon_and_4_5_point_towards_and_away_from_the_viewer():
    assert formula(edit(methane(), type="fragment", atom=0, key="1")["molecule"]) == "C2H6"
    for key, sign in (("4", 1), ("5", -1)):
        r = edit(methane(), type="fragment", atom=0, key=key)
        new_c = r["focus"]
        assert r["molecule"]["atoms"][new_c]["el"] == "C"
        assert sign * r["molecule"]["atoms"][new_c]["z"] > 0.5  # default view looks down -z


def test_group_key_on_an_internal_atom_adds_a_substituent():
    propane, _ = chem.from_smiles("CCC")
    r = edit(propane, type="fragment", atom=1, key="3")
    assert formula(r["molecule"]) == "C9H12"  # isopropylbenzene
    assert any("substituent" in n for n in r["notes"])
    sane(r["molecule"])


def test_label_key_on_terminal_carbon_replaces_it():
    m, c2 = ethane()
    r = edit(m, type="fragment", atom=c2, key="o")
    assert formula(r["molecule"]) == "CH4O"  # the CH3 *becomes* OH: methanol
    assert r["molecule"]["atoms"][r["focus"]]["el"] == "O"
    sane(r["molecule"])


def test_label_key_on_hydrogen_replaces_that_hydrogen():
    m = methane()
    r = edit(m, type="fragment", atom=2, key="o")
    assert formula(r["molecule"]) == "CH4O"
    assert r["map"][2] == -1
    sane(r["molecule"])


def test_nitro_group_keeps_charges_and_neutral_total():
    benzene, _ = chem.from_smiles("c1ccccc1")
    h = next(i for i, a in enumerate(benzene["atoms"]) if a["el"] == "H")
    m = edit(benzene, type="fragment", atom=h, key="N")["molecule"]
    assert formula(m) == "C6H5NO2"
    assert sum(a["charge"] for a in m["atoms"]) == 0
    sane(m)


def test_label_key_inside_a_ring_changes_the_element():
    ring, _ = chem.from_smiles("C1CCCCC1")
    r = edit(ring, type="fragment", atom=0, key="o")
    assert formula(r["molecule"]) == "C5H10O"
    sane(r["molecule"])
    cf3 = edit(ring, type="fragment", atom=0, key="F")
    assert formula(cf3["molecule"]) == "C7H11F3" and cf3["notes"]
    with pytest.raises(ChemError, match="cannot replace"):
        edit(ring, type="fragment", atom=0, key="f")


def test_label_key_on_isolated_atom():
    r = edit(methane(), type="fragment", atom=0, key="o")
    assert formula(r["molecule"]) == "H2O"


def test_chain_key_focus_moves_to_new_atom_and_vinyl_to_its_end():
    r = edit(methane(), type="fragment", atom=0, key="1")
    assert formula(r["molecule"]) == "C2H6" and r["molecule"]["atoms"][r["focus"]]["el"] == "C"
    r = edit(methane(), type="fragment", atom=0, key="8")
    m, f = r["molecule"], r["focus"]
    assert formula(m) == "C2H4"  # the carbon became a vinyl group: ethylene
    assert sum(1 for j in neighbors(m, f) if m["atoms"][j]["el"] == "H") == 2  # terminal =CH2


def test_fragments_without_hydrogens_in_skeleton_mode():
    bare = {"atoms": [{"el": "C", "x": 0, "y": 0, "z": 0, "charge": 0}], "bonds": []}
    assert formula(edit(bare, keep_h=False, type="fragment", atom=0, key="3")["molecule"]) == "C6"


def test_fragment_avoids_clashes_with_existing_atoms():
    # Attach tert-butyl next to an existing tert-butyl: no atoms on top of each other.
    m, _ = chem.from_smiles("CC(C)(C)C")
    h = next(i for i, a in enumerate(m["atoms"]) if a["el"] == "H")
    r = edit(m, type="fragment", atom=h, key="K")
    assert min_nonbonded(r["molecule"]) > 1.2


@pytest.mark.parametrize("frag", FRAGMENTS, ids=lambda f: f.key)
def test_every_fragment_on_a_methane_hydrogen_is_valid_and_optimizable(frag):
    r = edit(methane(), type="fragment", atom=1, key=frag.key)
    m = r["molecule"]
    sane(m)
    method = "MMFF94"
    try:
        result = chem.optimize(m, method)
    except ChemError:
        method = "UFF"
        result = chem.optimize(m, method)
    assert math.isfinite(result["energy"])
    assert r["focus"] is not None


@pytest.mark.parametrize("frag", FRAGMENTS, ids=lambda f: f.key)
def test_every_fragment_template_has_one_attachment(frag):
    t = template(frag.key)
    assert t.elements[t.dummy] == "H" and t.attach != t.dummy


def test_edit_errors():
    m = methane()
    with pytest.raises(ChemError, match="Unknown edit"):
        apply_edit(m, {"type": "explode"})
    with pytest.raises(ChemError, match="No fragment"):
        apply_edit(m, {"type": "fragment", "atom": 0, "key": "%"})
    with pytest.raises(ChemError, match="does not exist"):
        apply_edit(m, {"type": "grow", "atom": 9, "element": "C"})
    with pytest.raises(ChemError, match="view"):
        apply_edit(m, {"type": "grow", "atom": 0, "element": "C"}, view={"right": [0, 0, 0], "up": [0, 1, 0], "toCamera": [0, 0, 1]})
    with pytest.raises(ChemError, match="keep_h"):
        apply_edit(m, {"type": "grow", "atom": 0, "element": "C"}, keep_h="yes")


# ------------------------------------------------------------- ring fusion


def _cc_bond(mol, order=None):
    for b in mol["bonds"]:
        if mol["atoms"][b["a"]]["el"] == "C" and mol["atoms"][b["b"]]["el"] == "C" and (order is None or b["order"] == order):
            return b["a"], b["b"]
    raise AssertionError("no C-C bond")


@pytest.mark.parametrize(
    "smiles, size, aromatic, expected",
    [
        ("C1CCCCC1", 6, True, "C10H12"),  # tetralin
        ("c1ccccc1", 6, True, "C10H8"),  # naphthalene
        ("CC", 3, False, "C3H6"),  # cyclopropane
        ("CC", 6, False, "C6H12"),  # cyclohexane
        ("C1CCCC1", 5, False, "C8H14"),  # bicyclo[3.3.0]octane
    ],
)
def test_fuse_ring_onto_bond(smiles, size, aromatic, expected):
    mol, _ = chem.from_smiles(smiles)
    a, b = _cc_bond(mol)
    r = edit(mol, type="fuse", a=a, b=b, size=size, aromatic=aromatic)
    assert formula(r["molecule"]) == expected
    sane(r["molecule"])
    result = chem.optimize(r["molecule"], "MMFF94")
    assert math.isfinite(result["energy"])


def test_fused_benzene_is_aromatic_by_rdkit():
    from rdkit import Chem

    mol, _ = chem.from_smiles("c1ccccc1")
    a, b = _cc_bond(mol, order=1)
    r = edit(mol, type="fuse", a=a, b=b, size=6, aromatic=True)
    rd = chem.mol_from_json(r["molecule"])
    assert Chem.MolToSmiles(Chem.RemoveHs(rd)) == "c1ccc2ccccc2c1"


def test_fuse_errors():
    mol, _ = chem.from_smiles("C#C")
    a, b = _cc_bond(mol)
    with pytest.raises(ChemError, match="triple"):
        edit(mol, type="fuse", a=a, b=b, size=6, aromatic=True)
    m = methane()
    with pytest.raises(ChemError, match="heavy atoms"):
        edit(m, type="fuse", a=0, b=1, size=6)
    e, c2 = ethane()
    with pytest.raises(ChemError, match="six-membered"):
        edit(e, type="fuse", a=0, b=c2, size=5, aromatic=True)


def test_bond_key_table_matches_edits():
    from minimol.edit import BOND_KEYS

    for key, spec in BOND_KEYS.items():
        mol, _ = chem.from_smiles("CCCC")
        a, b = 1, 2  # the middle C-C bond
        if spec[0] == "order":
            r = edit(mol, type="bond_order", a=a, b=b, order=spec[1])
        else:
            r = edit(mol, type="fuse", a=a, b=b, size=spec[1], aromatic=spec[2])
        sane(r["molecule"])


def test_over_valent_atoms_are_reported_in_3d():
    m = {"atoms": [{"el": "C", "x": 0.0, "y": 0.0, "z": 0.0, "charge": 0}] + [
        {"el": "C", "x": 1.5 * x, "y": 1.5 * y, "z": 1.5 * z, "charge": 0}
        for x, y, z in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))
    ], "bonds": [{"a": 0, "b": k, "order": 1} for k in range(1, 5)]}
    r = edit(m, keep_h=False, type="grow", atom=0, element="C")
    assert any("Too many bonds on C1" in n for n in r["notes"])


# --------------------------------------- fusing next to existing rings (ChemDraw-like)


def _smiles(mol):
    from rdkit import Chem

    rd = chem.mol_from_json(mol, sanitize=False)
    for a in rd.GetAtoms():
        a.SetNoImplicit(a.GetAtomicNum() == 1)
    Chem.SanitizeMol(rd)  # also proves no atom is over its valence
    return Chem.MolToSmiles(Chem.RemoveHs(rd))


def _canon(smiles):
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def _heavy_degree(mol, i):
    return sum(1 for j in neighbors(mol, i) if mol["atoms"][j]["el"] != "H")


def _molecule(smiles, dim):
    from minimol import layout

    return chem.from_smiles(smiles)[0] if dim == 3 else layout.smiles_2d(smiles)["molecule"]


@pytest.mark.parametrize("dim", [3, 2])
def test_benzene_fused_onto_a_bond_between_two_rings_joins_both(dim):
    # The user's case: "a" on the bond linking two rings -> the new ring shares
    # atoms with both rings (biphenyl -> phenanthrene), not a ring stuck on the side.
    bp = _molecule("c1ccc(cc1)-c1ccccc1", dim)
    b = next(b for b in bp["bonds"] if _heavy_degree(bp, b["a"]) == 3 and _heavy_degree(bp, b["b"]) == 3)
    r = apply_edit(bp, {"type": "fuse", "a": b["a"], "b": b["b"], "size": 6, "aromatic": True}, dim=dim)
    assert _smiles(r["molecule"]) == _canon("c1ccc2c(c1)ccc1ccccc12")
    assert not r["notes"]


@pytest.mark.parametrize("dim", [3, 2])
def test_fusion_onto_a_bond_between_ring_and_heteroaromatic_ring(dim):
    m = _molecule("c1ccn(c1)-c1ccccc1", dim)
    b = next(
        b for b in m["bonds"]
        if {m["atoms"][b["a"]]["el"], m["atoms"][b["b"]]["el"]} == {"C", "N"}
        and _heavy_degree(m, b["a"]) == 3 and _heavy_degree(m, b["b"]) == 3
    )
    r = apply_edit(m, {"type": "fuse", "a": b["a"], "b": b["b"], "size": 6, "aromatic": True}, dim=dim)
    assert _smiles(r["molecule"]) == _canon("c1ccc2c(c1)ccc1cccn12")  # pyrrolo[1,2-a]quinoline


def test_fusion_onto_ring_ch2_bond_goes_through_the_ortho_carbon_in_3d():
    m = _molecule("CCc1ccccc1", 3)
    b = next(
        b for b in m["bonds"]
        if {m["atoms"][b["a"]]["el"], m["atoms"][b["b"]]["el"]} == {"C"}
        and sorted([_heavy_degree(m, b["a"]), _heavy_degree(m, b["b"])]) == [2, 3]
    )
    r = apply_edit(m, {"type": "fuse", "a": b["a"], "b": b["b"], "size": 6, "aromatic": True})
    assert _smiles(r["molecule"]) == _canon("Cc1cccc2ccccc12")  # 1-methylnaphthalene, nothing over-valent


def test_fusion_never_reuses_atoms_of_a_ring_that_contains_the_bond():
    for dim in (3, 2):
        m = _molecule("c1ccccc1", dim)
        b = m["bonds"][0]
        r = apply_edit(m, {"type": "fuse", "a": b["a"], "b": b["b"], "size": 6, "aromatic": True}, dim=dim)
        assert _smiles(r["molecule"]) == _canon("c1ccc2ccccc2c1")
