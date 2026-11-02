"""Geometry optimization with g-xTB through an external ``xtb`` program.

g-xTB (Grimme group, https://github.com/grimme-lab/g-xtb) is distributed as a
modified ``xtb`` binary that understands a ``--gxtb`` flag. It is licensed
under the GPL-3.0 and is *not* bundled with minimol: minimol only runs it as a
separate program when the user has installed it.

The program is looked up in this order:

1. the ``MINIMOL_GXTB`` environment variable (full path to the binary),
2. a ``gxtb`` program on ``PATH``,
3. an ``xtb`` program on ``PATH``.

A candidate is only used if its ``--help`` text lists the ``--gxtb`` flag, so a
regular xtb without g-xTB support is never mistaken for it.

Behaviour of the binary that this module relies on (checked against the
g-xTB 2.0.1 / xtb 6.7.1 build of 2026-05-14):

* ``xtb in.xyz --gxtb --opt`` writes the optimized structure to ``xtbopt.xyz``
  whose comment line reads ``energy: <Eh> gnorm: <Eh/bohr> ...``.
* It exits with status 0 both when the optimization converged
  (log: ``GEOMETRY OPTIMIZATION CONVERGED``) and when it ran out of cycles
  (log: ``FAILED TO CONVERGE``), so convergence is read from the log.
* Real errors give a non-zero exit status; the message is on log lines
  starting with ``-1-``, ``-2-``, ...
* ``$fix / atoms: 1,2 / $end`` in an xcontrol file (1-based indices) keeps those
  atoms exactly in place.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import tempfile
import threading

from .chem import ChemError, validate_json

ENV_VAR = "MINIMOL_GXTB"
MAX_ATOMIC_NUMBER = 103  # g-xTB covers H to Lr
OPT_LEVELS = ("crude", "sloppy", "loose", "normal", "tight", "verytight", "extreme")
DEFAULT_TIMEOUT = 1800  # seconds

_PERIODIC = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr"
).split()
_Z = {sym: i + 1 for i, sym in enumerate(_PERIODIC)}

_lock = threading.Lock()
_cache: dict[str, str | None] = {}


def _supports_gxtb(path: str) -> bool:
    try:
        result = subprocess.run(
            [path, "--help"],
            capture_output=True,
            text=True,
            timeout=30,
            stdin=subprocess.DEVNULL,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return re.search(r"^\s*--gxtb\b", result.stdout + result.stderr, re.M) is not None


def find_binary(refresh: bool = False) -> str | None:
    """Return the path of an xtb binary with g-xTB support, or ``None``."""
    with _lock:
        if not refresh and "path" in _cache:
            return _cache["path"]
        candidates = []
        env = os.environ.get(ENV_VAR)
        if env:
            candidates.append(os.path.expanduser(env))
        for name in ("gxtb", "xtb"):
            found = shutil.which(name)
            if found:
                candidates.append(found)
        path = None
        for candidate in candidates:
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK) and _supports_gxtb(candidate):
                path = candidate
                break
        _cache["path"] = path
        return path


def _electron_count(atoms: list[dict], charge: int) -> int:
    return sum(_Z[a["el"]] for a in atoms) - charge


def optimize(
    data: object,
    charge: object = None,
    uhf: object = None,
    fixed: object = (),
    level: object = "normal",
    max_cycles: object = None,
    timeout: float = DEFAULT_TIMEOUT,
    binary: str | None = None,
) -> dict:
    """Optimize ``data`` with g-xTB.

    ``charge`` defaults to the sum of the atoms' formal charges. ``uhf`` (the
    number of unpaired electrons) defaults to xtb's own choice: singlet for an
    even electron count, doublet for an odd one. Returns the same shape as
    :func:`minimol.chem.optimize`, with energies in Hartree.
    """
    validate_json(data)
    atoms = data["atoms"]
    n = len(atoms)
    if n == 0:
        raise ChemError("Nothing to optimize: the molecule is empty")
    for i, atom in enumerate(atoms):
        if _Z.get(atom["el"], MAX_ATOMIC_NUMBER + 1) > MAX_ATOMIC_NUMBER:
            raise ChemError(f"g-xTB covers H to Lr only; atom {i + 1} is {atom['el']}")

    if charge is None:
        charge = sum(a.get("charge", 0) for a in atoms)
    if isinstance(charge, bool) or not isinstance(charge, int) or abs(charge) > 50:
        raise ChemError("Total charge must be an integer")
    electrons = _electron_count(atoms, charge)
    if electrons < 0:
        raise ChemError(f"A total charge of {charge} leaves a negative number of electrons")
    if uhf is not None:
        if isinstance(uhf, bool) or not isinstance(uhf, int) or uhf < 0:
            raise ChemError("Unpaired electrons must be a non-negative integer")
        if uhf > electrons or (electrons - uhf) % 2:
            raise ChemError(
                f"{uhf} unpaired electron(s) is impossible with {electrons} electrons "
                f"(total charge {charge})"
            )
    if level not in OPT_LEVELS:
        raise ChemError(f"Optimization level must be one of {', '.join(OPT_LEVELS)}")
    if max_cycles is not None:
        if isinstance(max_cycles, bool) or not isinstance(max_cycles, int) or not 1 <= max_cycles <= 100000:
            raise ChemError("max_cycles must be between 1 and 100000")
    if not isinstance(fixed, (list, tuple)):
        raise ChemError("fixed must be a list of atom indices")
    fixed_atoms = sorted(set(fixed))
    if any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < n for i in fixed_atoms):
        raise ChemError("A fixed atom index is out of range")
    if len(fixed_atoms) == n:
        raise ChemError("All atoms are fixed; there is nothing to optimize")

    binary = binary or find_binary()
    if binary is None:
        raise ChemError(
            "g-xTB was not found. Install it from https://github.com/grimme-lab/g-xtb and put it "
            f"on PATH as 'gxtb', or set {ENV_VAR} to the binary's path."
        )

    with tempfile.TemporaryDirectory(prefix="minimol-gxtb-") as work:
        xyz_lines = [str(n), "minimol"]
        xyz_lines += [f"{a['el']} {a['x']:.10f} {a['y']:.10f} {a['z']:.10f}" for a in atoms]
        with open(os.path.join(work, "in.xyz"), "w") as fh:
            fh.write("\n".join(xyz_lines) + "\n")
        args = [binary, "in.xyz", "--gxtb", "--opt", level, "--norestart", "--chrg", str(charge)]
        if uhf is not None:
            args += ["--uhf", str(uhf)]
        if max_cycles is not None:
            args += ["--cycles", str(max_cycles)]
        if fixed_atoms:
            with open(os.path.join(work, "xcontrol"), "w") as fh:
                fh.write("$fix\n   atoms: " + ",".join(str(i + 1) for i in fixed_atoms) + "\n$end\n")
            args += ["--input", "xcontrol"]

        log_path = os.path.join(work, "xtb.log")
        try:
            with open(log_path, "w") as log:
                proc = subprocess.run(
                    args,
                    cwd=work,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,
                    timeout=timeout,
                )
        except subprocess.TimeoutExpired:
            raise ChemError(f"g-xTB did not finish within {timeout:g} seconds") from None
        except OSError as exc:
            raise ChemError(f"Could not run g-xTB ({binary}): {exc}") from None
        with open(log_path, errors="replace") as fh:
            log_text = fh.read()

        out_path = os.path.join(work, "xtbopt.xyz")
        if proc.returncode != 0 or not os.path.exists(out_path):
            raise ChemError("g-xTB failed: " + _error_message(log_text, proc.returncode))
        with open(out_path) as fh:
            optimized, energy, gnorm = _read_xtbopt(fh.read(), atoms)

    converged = "GEOMETRY OPTIMIZATION CONVERGED" in log_text
    notes = []
    if not converged:
        notes.append("g-xTB did not converge within the cycle limit; run Optimize again to continue.")
    initial = _initial_energy(log_text)
    molecule = {
        "atoms": [
            {**atom, "x": x, "y": y, "z": z}
            for atom, (x, y, z) in zip(atoms, optimized)
        ],
        "bonds": [dict(b) for b in data.get("bonds", [])],
    }
    return {
        "molecule": molecule,
        "energy": energy,
        "initial_energy": initial,
        "energy_unit": "Eh",
        "gradient_norm": gnorm,
        "converged": converged,
        "method": "g-xTB",
        "charge": charge,
        "uhf": uhf,
        "notes": notes,
    }


def _read_xtbopt(text: str, atoms: list[dict]) -> tuple[list[tuple[float, float, float]], float, float | None]:
    lines = text.splitlines()
    try:
        count = int(lines[0].split()[0])
    except (IndexError, ValueError):
        raise ChemError("g-xTB wrote an unreadable xtbopt.xyz") from None
    if count != len(atoms) or len(lines) < count + 2:
        raise ChemError("g-xTB returned a different number of atoms than it was given")
    match = re.search(r"energy:\s*(\S+)", lines[1])
    if not match:
        raise ChemError("g-xTB output has no energy")
    energy = float(match.group(1))
    gmatch = re.search(r"gnorm:\s*(\S+)", lines[1])
    gnorm = float(gmatch.group(1)) if gmatch else None
    coords = []
    for i, line in enumerate(lines[2 : count + 2]):
        parts = line.split()
        sym = parts[0][:1].upper() + parts[0][1:].lower()
        if sym != atoms[i]["el"]:
            raise ChemError(f"g-xTB changed the element of atom {i + 1} ({atoms[i]['el']} -> {sym})")
        xyz = tuple(float(v) for v in parts[1:4])
        if not all(math.isfinite(v) for v in xyz):
            raise ChemError("g-xTB returned non-finite coordinates")
        coords.append(tuple(round(v, 8) + 0.0 for v in xyz))  # + 0.0 turns -0.0 into 0.0
    if not math.isfinite(energy):
        raise ChemError("g-xTB returned a non-finite energy")
    return coords, energy, gnorm


def _initial_energy(log_text: str) -> float | None:
    # The first SCF of the run reports ":: total energy  <value> Eh ::".
    match = re.search(r"::\s*total energy\s+(-?\d+\.\d+)\s*Eh", log_text)
    return float(match.group(1)) if match else None


def _error_message(log_text: str, returncode: int) -> str:
    lines = [
        re.sub(r"^-\d+-\s*", "", line).strip()
        for line in log_text.splitlines()
        if re.match(r"^-\d+-\s", line)
    ]
    if lines:
        # xtb prints the innermost cause with the lowest number last.
        return "; ".join(reversed(lines))
    tail = [line.strip() for line in log_text.splitlines() if line.strip()][-3:]
    return f"exit status {returncode}" + (f" ({' / '.join(tail)})" if tail else "")
