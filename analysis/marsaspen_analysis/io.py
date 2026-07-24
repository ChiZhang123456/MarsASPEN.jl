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


def particle_history(data: dict[str, np.ndarray], particle_id: int) -> dict[str, np.ndarray]:
    """Return all event columns for one particle."""
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
    """Return collision rows, optionally restricted to chemical reactions."""
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
