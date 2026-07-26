"""Read the original processed GITM and AMPS MATLAB atmosphere files."""

from pathlib import Path

import numpy as np
from scipy.io import loadmat


def _case_names(ls: int) -> tuple[str, str]:
    cases = {
        0: ("aequ", "subL0"),
        90: ("aph", "subL0"),
        180: ("aequ", "subL180"),
        270: ("per", "subL270"),
    }
    if ls not in cases:
        raise ValueError("ls must be 0, 90, 180, or 270")
    return cases[ls]


def _cube(data: dict, field: str, shape: tuple[int, int, int]) -> np.ndarray:
    """Restore a MATLAB-vectorized field and remove duplicate longitude."""
    return np.asarray(data[field], dtype=float).reshape(shape, order="F")[:72]


def _read_gitm(path: Path) -> dict[str, np.ndarray]:
    data = loadmat(path, squeeze_me=True)
    shape = (73, 36, 50)
    longitude = _cube(data, "Longitude", shape)
    latitude = _cube(data, "Latitude", shape)
    altitude = _cube(data, "Altitude", shape)
    result = {
        "lon_deg": longitude[:, 0, 0],
        "lat_deg": latitude[0, :, 0],
        "alt_km": altitude[0, 0, :],
        "Tn": _cube(data, "temp", shape),
        "nCO2": _cube(data, "NCO2", shape) * 1e6,
        "nO": _cube(data, "NO", shape) * 1e6,
        "SZA_deg": _cube(data, "SZA", shape),
        "source_file": path.name,
    }
    floor = np.full(result["nCO2"].shape, 1e-300)
    result.update({"nO2": floor, "nN2": floor, "nCO": floor})
    return result


def _read_amps(path: Path) -> dict[str, np.ndarray]:
    data = loadmat(path, squeeze_me=True)
    shape = (73, 36, 167)
    longitude = _cube(data, "Longitude", shape)
    latitude = _cube(data, "Latitude", shape)
    altitude = _cube(data, "Altitude", shape)
    return {
        "lon_deg": longitude[:, 0, 0],
        "lat_deg": latitude[0, :, 0],
        "alt_km": altitude.mean(axis=(0, 1)),
        "nO_hot": _cube(data, "dens_oh", shape) * 1e6,
        "source_file": path.name,
    }


def _log_blend(low: np.ndarray, high: np.ndarray, weight: float) -> np.ndarray:
    return np.exp(
        (1.0 - weight) * np.log(np.maximum(low, 1e-300))
        + weight * np.log(np.maximum(high, 1e-300))
    )


def load_atmosphere_case(
    atmosphere_dir: str | Path, ls: int, f107: int
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Load one MSO GITM and AMPS case using the same rules as Julia."""
    if f107 not in (70, 130, 200):
        raise ValueError("f107 must be 70, 130, or 200")
    root = Path(atmosphere_dir)
    season, subsolar_label = _case_names(ls)

    def paths(activity: str) -> tuple[Path, Path]:
        return (
            root / "GITM" / f"gitm_{season}{activity}_{subsolar_label}_alt220.mat",
            root / "AMPS" / f"dsmc_{season}{activity}.mat",
        )

    low_paths, high_paths = paths("min"), paths("max")
    if f107 == 70:
        return _read_gitm(low_paths[0]), _read_amps(low_paths[1])
    if f107 == 200:
        return _read_gitm(high_paths[0]), _read_amps(high_paths[1])

    weight = (130.0 - 70.0) / (200.0 - 70.0)
    gl, al = _read_gitm(low_paths[0]), _read_amps(low_paths[1])
    gh, ah = _read_gitm(high_paths[0]), _read_amps(high_paths[1])
    gitm = {
        "lon_deg": gl["lon_deg"],
        "lat_deg": gl["lat_deg"],
        "alt_km": gl["alt_km"],
        "Tn": (1.0 - weight) * gl["Tn"] + weight * gh["Tn"],
        "nCO2": _log_blend(gl["nCO2"], gh["nCO2"], weight),
        "nO": _log_blend(gl["nO"], gh["nO"], weight),
        "SZA_deg": gl["SZA_deg"],
        "source_file": f"{gl['source_file']} + {gh['source_file']}",
    }
    floor = np.full(gitm["nCO2"].shape, 1e-300)
    gitm.update({"nO2": floor, "nN2": floor, "nCO": floor})
    amps = {
        "lon_deg": al["lon_deg"],
        "lat_deg": al["lat_deg"],
        "alt_km": al["alt_km"],
        "nO_hot": _log_blend(al["nO_hot"], ah["nO_hot"], weight),
        "source_file": f"{al['source_file']} + {ah['source_file']}",
    }
    return gitm, amps
