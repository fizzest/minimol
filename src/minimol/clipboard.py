"""Copy and paste structures through the system clipboard (e.g. with ChemDraw).

Paste looks at every format on the clipboard and recognises the data by its
content, not only by the format's name, because the names differ between
operating systems:

* CDXML (ChemDraw's XML format), read with RDKit;
* binary CDX (ChemDraw's native format, starts with "VjCD0100") when the
  installed RDKit can read it; the RDKit 2026.09 release cannot (its Python
  reader is fixed to CDXML), and you are told so;
* MDL CT ("MDLCT", used by ChemDraw and ISIS/Draw on Windows): a molfile in
  which every line is preceded by one byte holding its length;
* molfile text, SMILES and InChI text.

Copy puts several formats on the clipboard: SMILES as plain text (paste it
into ChemDraw with Paste Special > SMILES), a molfile, CDXML, and an MDL CT
version of the molfile under the Windows format name "MDLCT".
"""

from __future__ import annotations

import os
import tempfile
from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import rdChemDraw

from . import layout
from .chem import ChemError, _embed_3d, mol_from_json, mol_to_json, validate_json

CDX_MAGIC = b"VjCD0100"
MDLCT_WINDOWS_FORMAT = 'application/x-qt-windows-mime;value="MDLCT"'
COPY_FORMATS = ("text/plain", "chemical/x-mdl-molfile", "chemical/x-cdxml", MDLCT_WINDOWS_FORMAT)


# --------------------------------------------------------------------------
# MDL CT encoding
# --------------------------------------------------------------------------


def mdlct_encode(molblock: str) -> bytes:
    """Encode a molfile as MDL CT: each line preceded by a byte with its length."""
    out = bytearray()
    for line in molblock.splitlines():
        raw = line.encode("latin-1", errors="replace")[:255]
        out.append(len(raw))
        out += raw
    return bytes(out)


def mdlct_decode(data: bytes) -> str | None:
    """Decode MDL CT bytes to molfile text, or None if the data is not MDL CT."""
    data = data.rstrip(b"\x00")
    lines = []
    i = 0
    while i < len(data):
        n = data[i]
        i += 1
        if i + n > len(data):
            return None
        chunk = data[i : i + n]
        if b"\n" in chunk or b"\r" in chunk:
            return None
        lines.append(chunk.decode("latin-1"))
        i += n
    text = "\n".join(lines)
    return text + "\n" if "M  END" in text else None


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------


@lru_cache(maxsize=1)
def cdx_supported() -> bool:
    """Whether this RDKit can read binary CDX (needs the ``format`` argument, added after 2026.09)."""
    try:
        rdChemDraw.MolsFromChemDrawFile(os.path.join(tempfile.gettempdir(), "minimol-probe.cdx"), format=rdChemDraw.CDXFormat.CDX)
    except TypeError:  # Boost.Python's ArgumentError: no "format" argument in this RDKit
        return False
    except Exception:  # the probe file does not exist; the argument itself was accepted
        return True
    return True


def _combine(mols: list[Chem.Mol]) -> Chem.Mol | None:
    mols = [m for m in mols if m is not None and m.GetNumAtoms()]
    if not mols:
        return None
    out = mols[0]
    for m in mols[1:]:
        out = Chem.CombineMols(out, m)
    return out


def _from_cdxml(text: str) -> Chem.Mol | None:
    try:
        return _combine(list(Chem.MolsFromCDXML(text, removeHs=False)))
    except Exception:
        return None


def _from_cdx(data: bytes) -> Chem.Mol | None:
    fd, path = tempfile.mkstemp(suffix=".cdx")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        mols = rdChemDraw.MolsFromChemDrawFile(path, format=rdChemDraw.CDXFormat.CDX)
        return _combine(list(mols))
    except Exception:
        return None
    finally:
        os.unlink(path)


def _from_molblock(text: str) -> Chem.Mol | None:
    try:
        return Chem.MolFromMolBlock(text, removeHs=False)
    except Exception:
        return None


def _from_text(raw: str) -> tuple[Chem.Mol | None, str]:
    raw = raw.strip("\x00")
    if "M  END" in raw:
        # Not stripped: a molfile's first line (the title) may be blank.
        return _from_molblock(raw.rstrip() + "\n"), "molfile text"
    text = raw.strip()
    if not text:
        return None, ""
    if "<CDXML" in text[:2000]:
        return _from_cdxml(text), "CDXML text"
    if text.startswith("InChI="):
        try:
            return Chem.MolFromInchi(text.split()[0]), "InChI"
        except Exception:
            return None, "InChI"
    first = text.split()[0] if len(text.split()) else ""
    if len(text.splitlines()) <= 2 and first:
        mol = Chem.MolFromSmiles(first)
        if mol is not None:
            return mol, "SMILES"
    return None, ""


def read_structure(formats: dict[str, bytes]) -> tuple[Chem.Mol, str]:
    """Find a structure in the clipboard formats. Returns (molecule, description of the source)."""
    saw_cdx = False
    # 1. Rich chemical formats, recognised by content.
    for name, data in formats.items():
        if not data:
            continue
        head = data[:4000]
        if b"<CDXML" in head:
            mol = _from_cdxml(data.decode("utf-8", errors="replace"))
            if mol is not None:
                return mol, f"CDXML ({name})"
        if data.startswith(CDX_MAGIC):
            saw_cdx = True
            if cdx_supported():
                mol = _from_cdx(data)
                if mol is not None:
                    return mol, f"ChemDraw CDX ({name})"
    for name, data in formats.items():
        if not data:
            continue
        if "mdlct" in name.lower() or (data[0] < 0x50 and b"M  END" in data and b"\n" not in data[:200]):
            text = mdlct_decode(data)
            if text:
                mol = _from_molblock(text)
                if mol is not None:
                    return mol, f"MDL CT ({name})"
    for name, data in formats.items():
        if data and ("molfile" in name.lower() or "mdl" in name.lower()) and b"M  END" in data:
            mol = _from_molblock(data.decode("latin-1"))
            if mol is not None:
                return mol, f"molfile ({name})"
    # 2. Plain text: molfile, CDXML, InChI or SMILES.
    for name, data in formats.items():
        if data and (name.startswith("text/plain") or name.lower() in ("text", "utf8_string", "string")):
            mol, kind = _from_text(data.decode("utf-8", errors="replace"))
            if mol is not None:
                return mol, kind
    listed = ", ".join(sorted(formats)) or "none"
    if saw_cdx:
        raise ChemError(
            "The clipboard holds ChemDraw's binary CDX format, which the installed RDKit cannot read yet. "
            "In ChemDraw, copy the structure as SMILES or MOL text instead (Edit > Copy As), or enable "
            f"MDL CT clipboard output. Formats found: {listed}"
        )
    raise ChemError(f"Nothing on the clipboard could be read as a structure. Formats found: {listed}")


def paste(formats: dict[str, bytes], dim: int) -> dict:
    """Turn clipboard contents into a molecule for the editor (2D drawing or 3D structure)."""
    mol, source = read_structure(formats)
    notes = []
    if dim == 2:
        molecule = layout.from_rdkit_2d(mol, notes)
    else:
        work = Chem.RWMol(layout.strip_hydrogens(mol))
        for a in work.GetAtoms():
            a.SetNoImplicit(a.GetAtomicNum() == 1)
        try:
            Chem.SanitizeMol(work)
        except Exception as exc:
            raise ChemError(f"The pasted structure is not valid: {exc}") from None
        work.RemoveAllConformers()
        molecule = mol_to_json(_embed_3d(work), notes)
    return {"molecule": molecule, "notes": notes, "source": source}


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------


def _subset(data: dict, atoms: list[int] | None) -> Chem.RWMol:
    mol = mol_from_json(data, sanitize=False)
    if atoms:
        keep = set(atoms)
        if not all(isinstance(i, int) and 0 <= i < mol.GetNumAtoms() for i in keep):
            raise ChemError("A selected atom does not exist")
        for idx in sorted((i for i in range(mol.GetNumAtoms()) if i not in keep), reverse=True):
            mol.RemoveAtom(idx)
    return mol


def payload(data: object, atoms: object = None, dim: int = 3) -> dict:
    """Clipboard representations of the molecule (or of the selected atoms).

    Structures always go out as 2D drawings with implicit hydrogens, which is
    what drawing programs expect: a 2D molecule keeps its layout, a 3D one is
    laid out flat.
    """
    validate_json(data)
    if atoms is not None and not isinstance(atoms, list):
        raise ChemError("atoms must be a list")
    mol = layout.strip_hydrogens(_subset(data, atoms or None))
    if mol.GetNumAtoms() == 0:
        raise ChemError("Nothing to copy")
    for a in mol.GetAtoms():
        a.SetNoImplicit(a.GetAtomicNum() == 1)
    mol.UpdatePropertyCache(strict=False)
    if dim == 2:
        mol.GetConformer().Set3D(False)  # already a drawing at the standard bond length
    else:
        layout.depict(mol)
    try:
        Chem.SanitizeMol(mol)
        smiles = Chem.MolToSmiles(mol)
    except Exception:
        smiles = Chem.MolToSmiles(mol, canonical=False)
    molblock = Chem.MolToMolBlock(mol)
    try:
        cdxml = Chem.MolToCDXMLBlock(mol)
    except Exception:
        cdxml = ""
    return {"smiles": smiles, "molblock": molblock, "cdxml": cdxml}


def copy_formats(p: dict) -> dict[str, bytes]:
    out = {
        "text/plain": p["smiles"].encode(),
        "chemical/x-mdl-molfile": p["molblock"].encode(),
        MDLCT_WINDOWS_FORMAT: mdlct_encode(p["molblock"]),
    }
    if p["cdxml"]:
        out["chemical/x-cdxml"] = p["cdxml"].encode()
    return out
