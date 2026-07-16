"""Partial charge assignment for electrostatic grids (PACMOF DDEC)."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional, Union

import numpy as np

logger = logging.getLogger(__name__)

CHARGE_INFO_KEY = "_atom_site_charge"
CHARGE_INFO_KEY_ALT = "_atom_site_charges"

# Fallback when pip install omits the trained pickle (known pacmof packaging issue).
DEFAULT_PACMOF_REPO_MODEL = (
    Path.home() / "repos" / "pacmof" / "pacmof" / "data" / "Model_RF_DDEC.pkl"
)


def resolve_pacmof_model_path(
    path_to_pickle_obj: Optional[Union[str, Path]] = None,
) -> Path:
    """Locate ``Model_RF_DDEC.pkl`` for PACMOF charge assignment.

    Search order:
    1. Explicit ``path_to_pickle_obj``
    2. Environment variable ``PACMOF_MODEL_PATH``
    3. ``~/repos/pacmof/pacmof/data/Model_RF_DDEC.pkl``
    4. ``pacmof`` package ``data/Model_RF_DDEC.pkl``
    """
    if path_to_pickle_obj is not None:
        p = Path(path_to_pickle_obj).expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(f"PACMOF model not found: {p}")
        return p

    candidates: list[Path] = []
    env = os.environ.get("PACMOF_MODEL_PATH")
    if env:
        candidates.append(Path(env).expanduser())
    candidates.append(DEFAULT_PACMOF_REPO_MODEL)
    try:
        import pacmof as pac

        candidates.append(Path(pac.__file__).resolve().parent / "data" / "Model_RF_DDEC.pkl")
    except ImportError:
        pass

    for p in candidates:
        if p.is_file():
            logger.debug("Using PACMOF model: %s", p)
            return p.resolve()

    searched = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        "PACMOF Model_RF_DDEC.pkl not found. Tried:\n  "
        f"{searched}\n"
        "Set PACMOF_MODEL_PATH, pass path_to_pickle_obj, or run pacmof/data/model_generator.py."
    )


def normalize_charge_info(atoms) -> np.ndarray:
    """Return per-atom charges in e, ensuring ``atoms.info['_atom_site_charge']``."""
    if CHARGE_INFO_KEY in atoms.info:
        charges = np.asarray(atoms.info[CHARGE_INFO_KEY], dtype=np.float64)
    elif CHARGE_INFO_KEY_ALT in atoms.info:
        charges = np.asarray(atoms.info[CHARGE_INFO_KEY_ALT], dtype=np.float64)
        atoms.info[CHARGE_INFO_KEY] = charges.tolist()
    else:
        raise KeyError(
            f"No charges on atoms. Expected info['{CHARGE_INFO_KEY}'] "
            f"(PACMOF sets this after get_charges_single_serial)."
        )
    if len(charges) != len(atoms):
        raise ValueError(f"Charge length {len(charges)} != natoms {len(atoms)}")
    return charges


def get_pacmof_charges_single(
    path_to_cif: Union[str, Path],
    *,
    path_to_pickle_obj: Optional[Union[str, Path]] = None,
    create_cif: bool = False,
) -> "ase.Atoms":
    """Assign DDEC-style partial charges with PACMOF.

    Returns ASE ``Atoms`` with neutralized charges in
    ``atoms.info['_atom_site_charge']`` (elementary charge, Kelvin/RASPA convention).

    The sklearn model path is resolved via :func:`resolve_pacmof_model_path` unless
    ``path_to_pickle_obj`` is given. PACMOF's broken ``use_default_model=True`` install
    path is avoided — we always pass an explicit pickle path.
    """
    import pacmof as pac

    path_to_cif = str(Path(path_to_cif).resolve())
    model_path = resolve_pacmof_model_path(path_to_pickle_obj)

    atoms = pac.get_charges_single_serial(
        path_to_cif,
        create_cif=create_cif,
        use_default_model=False,
        path_to_pickle_obj=str(model_path),
    )
    charges = normalize_charge_info(atoms)
    total = float(charges.sum())
    if abs(total) > 1e-4:
        logger.warning("PACMOF charges sum to %.2e e after neutrality adjustment.", total)
    else:
        logger.info(
            "PACMOF charges assigned: %d atoms, |sum q|=%.2e (model=%s)",
            len(atoms),
            abs(total),
            model_path.name,
        )
    return atoms


def atoms_with_pacmof_charges(path_to_cif: Union[str, Path], **kwargs):
    """Alias for :func:`get_pacmof_charges_single`."""
    return get_pacmof_charges_single(path_to_cif, **kwargs)
