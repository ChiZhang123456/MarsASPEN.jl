"""Unit tests for MAT v7.3 loading and event selectors."""

from pathlib import Path

import h5py
import numpy as np

from marsaspen_analysis import (
    load_history_mat,
    mat_string,
    particle_history,
    reaction_events,
)


def test_hdf5_history_reader(tmp_path: Path) -> None:
    filename = tmp_path / "history.mat"
    with h5py.File(filename, "w") as handle:
        handle["particle_id"] = np.array([1, 1, 2])
        handle["event_type"] = np.array([0, 2, 2], dtype=np.uint8)
        handle["reaction_code"] = np.array([0, 2, 4], dtype=np.uint8)
        handle["energy_ev"] = np.array([800.0, 780.0, 700.0])

    data = load_history_mat(filename)
    assert particle_history(data, 1)["energy_ev"].tolist() == [800.0, 780.0]
    assert reaction_events(data)["reaction_code"].tolist() == [2]
    assert reaction_events(data, include_elastic=True)["reaction_code"].tolist() == [2, 4]


def test_mat_string_decodes_uint16_characters() -> None:
    assert mat_string(np.array([72, 112, 108, 117, 115], dtype=np.uint16)) == "Hplus"
