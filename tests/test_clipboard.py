import pytest
from rdkit import Chem

from minimol import chem, clipboard, layout
from minimol.chem import ChemError

ETHANOL_MOL = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))


def smiles(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def test_mdlct_round_trip_and_format():
    data = clipboard.mdlct_encode(ETHANOL_MOL)
    # Each line is one length byte followed by the line, with no newlines.
    first_len = data[0]
    assert data[1 : 1 + first_len].decode() == ETHANOL_MOL.splitlines()[0]
    assert b"\n" not in data
    decoded = clipboard.mdlct_decode(data)
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(decoded)) == "CCO"
    assert clipboard.mdlct_decode(b"\x05ab") is None  # truncated
    assert clipboard.mdlct_decode(b"plain text\n") is None


@pytest.mark.parametrize(
    "formats, expected, source",
    [
        ({"text/plain": b"c1ccccc1O"}, "Oc1ccccc1", "SMILES"),
        ({"text/plain": b"InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"}, "CCO", "InChI"),
        ({"text/plain": ETHANOL_MOL.encode()}, "CCO", "molfile text"),
        ({"chemical/x-mdl-molfile": ETHANOL_MOL.encode()}, "CCO", "molfile"),
        ({'application/x-qt-windows-mime;value="MDLCT"': clipboard.mdlct_encode(ETHANOL_MOL)}, "CCO", "MDL CT"),
        ({"chemical/x-cdxml": Chem.MolToCDXMLBlock(Chem.MolFromSmiles("CCN")).encode()}, "CCN", "CDXML"),
        # Unnamed MDL CT data is recognised by content.
        ({"application/x-unknown": clipboard.mdlct_encode(ETHANOL_MOL)}, "CCO", "MDL CT"),
    ],
)
def test_read_structure_formats(formats, expected, source):
    mol, found = clipboard.read_structure(formats)
    assert smiles(mol) == expected
    assert found.startswith(source)


def test_rich_formats_win_over_text():
    formats = {
        "text/plain": b"C",
        "chemical/x-cdxml": Chem.MolToCDXMLBlock(Chem.MolFromSmiles("CCCl")).encode(),
    }
    mol, found = clipboard.read_structure(formats)
    assert smiles(mol) == "CCCl" and found.startswith("CDXML")


def test_unreadable_clipboard_lists_formats_and_explains_cdx():
    with pytest.raises(ChemError, match="Formats found: image/png, text/plain"):
        clipboard.read_structure({"text/plain": b"hello world, not chemistry", "image/png": b"\x89PNG"})
    if not clipboard.cdx_supported():
        with pytest.raises(ChemError, match="binary CDX"):
            clipboard.read_structure({"ChemDraw Interchange Format": clipboard.CDX_MAGIC + b"\x00" * 40})


def test_paste_as_2d_and_3d():
    r2 = clipboard.paste({"text/plain": b"CC(=O)O"}, dim=2)
    assert all(a["z"] == 0 for a in r2["molecule"]["atoms"]) and len(r2["molecule"]["atoms"]) == 4
    r3 = clipboard.paste({"text/plain": b"CC(=O)O"}, dim=3)
    assert len(r3["molecule"]["atoms"]) == 8 and r3["source"] == "SMILES"


def test_payload_whole_molecule_and_selection():
    ethanol, _ = chem.from_smiles("CCO")
    p = clipboard.payload(ethanol, None, dim=3)
    assert p["smiles"] == "CCO"
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(p["molblock"])) == "CCO"
    assert Chem.MolToSmiles(Chem.MolsFromCDXML(p["cdxml"])[0]) == "CCO"
    # The molfile is a flat drawing even for a 3D molecule.
    conf = Chem.MolFromMolBlock(p["molblock"]).GetConformer()
    assert all(abs(conf.GetAtomPosition(i).z) < 1e-6 for i in range(3))
    # Only the selected heavy atoms (C-O): methanol.
    o = next(i for i, a in enumerate(ethanol["atoms"]) if a["el"] == "O")
    c = 1
    assert clipboard.payload(ethanol, [c, o], dim=3)["smiles"] == "CO"
    two = layout.smiles_2d("c1ccccc1")["molecule"]
    assert clipboard.payload(two, None, dim=2)["smiles"] == "c1ccccc1"


def test_copy_formats_include_text_molfile_cdxml_and_mdlct():
    ethanol, _ = chem.from_smiles("CCO")
    formats = clipboard.copy_formats(clipboard.payload(ethanol, None, dim=3))
    assert formats["text/plain"] == b"CCO"
    assert set(formats) == set(clipboard.COPY_FORMATS)
    # What we copy, we can paste back.
    mol, _ = clipboard.read_structure(formats)
    assert smiles(mol) == "CCO"
