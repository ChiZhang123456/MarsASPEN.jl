"""Readers and selectors for MarsASPEN detailed MAT trajectory output."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat


def load_history_mat(filename: str | Path) -> dict[str, np.ndarray]:
    """Load a MarsASPEN detailed MAT file, including MAT v7.3 output."""
    path = Path(filename)
    try:
        return {
            key: np.asarray(value).squeeze()
            for key, value in loadmat(path, squeeze_me=True).items()
            if not key.startswith("__")
        }
    except (NotImplementedError, ValueError):
        with h5py.File(path, "r") as handle:
            return {
                key: np.asarray(value).squeeze()
                for key, value in handle.items()
                if isinstance(value, h5py.Dataset)
            }


def mat_string(value: np.ndarray | str) -> str:
    """Decode a MATLAB or Julia MAT string into ordinary Python text."""
    array = np.asarray(value).squeeze()
    if array.dtype.kind in "ui" and array.ndim == 1:
        return "".join(chr(int(code)) for code in array if int(code) != 0)
    if array.dtype.kind in "SU":
        return "".join(str(item) for item in np.ravel(array))
    return str(array)


def particle_history(data: dict[str, np.ndarray], particle_id: int) -> dict[str, np.ndarray]:
    """Return aligned one-dimensional event columns for one particle."""
    mask = np.asarray(data["particle_id"], dtype=int) == int(particle_id)
    n_events = mask.size
    return {
        key: np.asarray(value)[mask]
        for key, value in data.items()
        if np.asarray(value).ndim == 1 and np.asarray(value).size == n_events
    }


def reaction_events(
    data: dict[str, np.ndarray],
    particle_id: int | None = None,
    include_elastic: bool = False,
) -> dict[str, np.ndarray]:
    """Return collision rows, optionally restricted to chemical reactions.

    Reaction codes 1 to 3 are state change, ionization, and Ly-alpha.
    Elastic collisions use code 4 and are excluded unless requested.
    """
    event_type = np.asarray(data["event_type"], dtype=int)
    reaction = np.asarray(data["reaction_code"], dtype=int)
    mask = event_type == 2
    if not include_elastic:
        mask &= (reaction >= 1) & (reaction <= 3)
    if particle_id is not None:
        mask &= np.asarray(data["particle_id"], dtype=int) == int(particle_id)
    n_events = mask.size
    return {
        key: np.asarray(value)[mask]
        for key, value in data.items()
        if np.asarray(value).ndim == 1 and np.asarray(value).size == n_events
    }
