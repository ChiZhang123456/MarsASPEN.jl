"""Plot altitude, energy, charge state, speed, and reactions for one particle.

Julia saves both propagation rows and collision rows. ``event_type == 2``
identifies collisions, while ``reaction_code`` distinguishes state change,
target ionization, Ly-alpha production, and elastic scattering. Identical
reaction colors are used in the altitude and energy panels so each event can
be followed between physical quantities.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis.io import load_history_mat  # noqa: E402

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"

# State change means H+ charge exchange or neutral-H stripping, depending on
# the projectile charge before the collision.
REACTION_LABELS = {
    1: "State change",
    2: "Ionization",
    3: "Ly-alpha",
    4: "Elastic",
}
REACTION_COLORS = {
    1: "#d62728",
    2: "#ff7f0e",
    3: "#2ca02c",
    4: "#4d4d4d",
}


def vector(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    """Return one MAT column as a one-dimensional NumPy array."""
    return np.atleast_1d(np.asarray(data[key]).squeeze())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    # The loader supports classic MAT and HDF5-backed MAT v7.3. Event columns
    # remain aligned row by row.
    data = load_history_mat(args.mat_file)
    time_s = vector(data, "time_s")
    altitude = vector(data, "altitude_km")
    energy = vector(data, "energy_ev")
    charge = vector(data, "charge_state").astype(int)
    event_type = vector(data, "event_type").astype(int)
    reaction = vector(data, "reaction_code").astype(int)
    vx = vector(data, "vx_m_s")
    vy = vector(data, "vy_m_s")
    vz = vector(data, "vz_m_s")
    # Convert the SI velocity magnitude to km/s only for display.
    speed = np.sqrt(vx**2 + vy**2 + vz**2) / 1000

    initial_charge = int(vector(data, "config_initial_charge_state")[0])
    species = r"H$^+$" if initial_charge == 1 else "H ENA"
    # Propagation rows draw the tracks. Only collision rows receive markers.
    collision = event_type == 2

    fig, axes = plt.subplots(
        4, 1, figsize=(10, 12), sharex=True, constrained_layout=True
    )
    axes[0].plot(time_s, altitude, color="black", lw=1.2)
    axes[1].plot(time_s, energy, color="#1f77b4", lw=1.2)
    axes[2].step(time_s, charge, where="post", color="#9467bd", lw=1.2)
    axes[3].plot(time_s, speed, color="#8c564b", lw=1.2)

    # Elastic collisions are often numerous, so their markers are smaller.
    for code, label in REACTION_LABELS.items():
        mask = collision & (reaction == code)
        if not np.any(mask):
            continue
        size = 15 if code == 4 else 38
        axes[0].scatter(
            time_s[mask], altitude[mask], s=size,
            color=REACTION_COLORS[code], label=label, zorder=3,
        )
        axes[1].scatter(
            time_s[mask], energy[mask], s=size,
            color=REACTION_COLORS[code], zorder=3,
        )

    axes[0].set_ylabel("Altitude (km)")
    axes[1].set_ylabel("Energy (eV)")
    axes[2].set_ylabel("Charge state")
    axes[2].set_yticks([0, 1], ["H ENA", "H+"])
    axes[3].set_ylabel("Speed (km/s)")
    axes[3].set_xlabel("Time (s)")
    axes[0].legend(ncol=2, fontsize=9)
    for ax in axes:
        ax.grid(True, color="0.85", lw=0.7)
    fig.suptitle(f"Single {species}, 600 km, 400 km/s")

    output = args.output or args.mat_file.with_suffix(".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()
