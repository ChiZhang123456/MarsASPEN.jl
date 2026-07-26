"""Create a publication-style figure for one simulated particle trajectory.

Julia saves both propagation rows and collision rows. ``event_type == 2``
identifies collisions, while ``reaction_code`` distinguishes state change,
target ionization, Ly-alpha production, and elastic scattering. Identical
reaction colors are used in the altitude and energy panels so each event can
be followed between physical quantities.

The displayed altitude interval is configurable. The five panels provide
complementary evidence: reaction altitude, kinetic-energy evolution, charge
state, particle speed or reaction timeline, and collision scattering angle
over the same selected time interval.

The output argument specifies the PNG file saved for the example.
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
mpl.rcParams["font.size"] = 8
mpl.rcParams["axes.linewidth"] = 0.8
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["legend.frameon"] = False

# State change means H+ charge exchange or neutral-H stripping, depending on
# the projectile charge before the collision.
REACTION_LABELS = {
    1: "State change",
    2: "Ionization",
    3: "Ly-alpha",
    4: "Elastic",
}
REACTION_COLORS = {
    1: "#C44E52",
    2: "#DD8452",
    3: "#55A868",
    4: "#6B6B6B",
}
REACTION_TIMELINE_LEVELS = {
    4: 0,  # Elastic
    2: 1,  # Target ionization
    3: 2,  # Ly-alpha emission
    1: 3,  # Projectile charge-state change
}


def vector(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    """Return one MAT column as a one-dimensional NumPy array."""
    return np.atleast_1d(np.asarray(data[key]).squeeze())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--altitude-min", type=float, default=80.0)
    parser.add_argument("--altitude-max", type=float, default=400.0)
    parser.add_argument(
        "--panel-d",
        choices=("speed", "reactions"),
        default="speed",
        help="Quantity shown in panel d.",
    )
    args = parser.parse_args()
    if args.altitude_min >= args.altitude_max:
        raise ValueError("--altitude-min must be smaller than --altitude-max.")

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
    vx_before = vector(data, "vx_before_m_s")
    vy_before = vector(data, "vy_before_m_s")
    vz_before = vector(data, "vz_before_m_s")
    # Convert the SI velocity magnitude to km/s only for display.
    speed = np.sqrt(vx**2 + vy**2 + vz**2) / 1000
    speed_before = np.sqrt(
        vx_before**2 + vy_before**2 + vz_before**2
    )
    speed_after = np.sqrt(vx**2 + vy**2 + vz**2)
    denominator = speed_before * speed_after
    cosine_angle = np.divide(
        vx_before * vx + vy_before * vy + vz_before * vz,
        denominator,
        out=np.ones_like(denominator, dtype=float),
        where=denominator > 0,
    )
    scattering_angle_deg = np.degrees(
        np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    )

    initial_charge = int(vector(data, "config_initial_charge_state")[0])
    species = r"H$^+$" if initial_charge == 1 else "H ENA"

    # Restrict every panel to the same atmospheric part of the trajectory.
    # Applying one common mask preserves alignment among time, position,
    # velocity, charge, and reaction columns.
    altitude_window = (
        (altitude >= args.altitude_min)
        & (altitude <= args.altitude_max)
    )
    if not np.any(altitude_window):
        raise ValueError(
            "This trajectory contains no samples from "
            f"{args.altitude_min:g} to {args.altitude_max:g} km."
        )
    time_s = time_s[altitude_window]
    altitude = altitude[altitude_window]
    energy = energy[altitude_window]
    charge = charge[altitude_window]
    event_type = event_type[altitude_window]
    reaction = reaction[altitude_window]
    speed = speed[altitude_window]
    scattering_angle_deg = scattering_angle_deg[altitude_window]

    # Propagation rows draw the tracks. Only collision rows receive markers.
    collision = event_type == 2

    fig = plt.figure(figsize=(7.2, 7.0), constrained_layout=True)
    grid = fig.add_gridspec(
        3, 2, height_ratios=(1.0, 1.0, 0.92)
    )
    axes = np.array(
        [
            fig.add_subplot(grid[0, 0]),
            fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, 0]),
            fig.add_subplot(grid[1, 1]),
            fig.add_subplot(grid[2, :]),
        ],
        dtype=object,
    )
    for axis in axes[1:]:
        axis.sharex(axes[0])
    axes[0].plot(time_s, altitude, color="#222222", lw=1.1)
    axes[1].plot(time_s, energy, color="#4C72B0", lw=1.1)
    axes[2].step(time_s, charge, where="post", color="#8172B2", lw=1.1)
    if args.panel_d == "speed":
        axes[3].plot(time_s, speed, color="#937860", lw=1.1)

    # Elastic collisions are often numerous, so their markers are smaller.
    for code, label in REACTION_LABELS.items():
        mask = collision & (reaction == code)
        if not np.any(mask):
            continue
        size = 9 if code == 4 else 24
        axes[0].scatter(
            time_s[mask], altitude[mask], s=size,
            color=REACTION_COLORS[code], label=label, zorder=3,
            edgecolors="none",
        )
        axes[1].scatter(
            time_s[mask], energy[mask], s=size,
            color=REACTION_COLORS[code], zorder=3, edgecolors="none",
        )
        axes[4].scatter(
            time_s[mask], scattering_angle_deg[mask], s=size,
            color=REACTION_COLORS[code], zorder=3, edgecolors="none",
        )

    axes[0].set_ylabel("Altitude (km)")
    axes[0].set_ylim(args.altitude_min, args.altitude_max)
    axes[1].set_ylabel("Energy (eV)")
    axes[2].set_ylabel("Charge state")
    axes[2].set_yticks([0, 1], ["H ENA", "H+"])
    if args.panel_d == "speed":
        axes[3].set_ylabel("Speed (km/s)")
    else:
        counts = {}
        for code, level in REACTION_TIMELINE_LEVELS.items():
            mask = collision & (reaction == code)
            counts[code] = int(np.count_nonzero(mask))
            if np.any(mask):
                axes[3].scatter(
                    time_s[mask],
                    np.full(np.count_nonzero(mask), level),
                    s=18 if code == 4 else 28,
                    color=REACTION_COLORS[code],
                    edgecolors="none",
                    zorder=3,
                )
        axes[3].set_ylabel("Reaction type")
        axes[3].set_yticks(
            [0, 1, 2, 3],
            ["Elastic", "Ionization", "Ly-alpha", "State change"],
        )
        axes[3].set_ylim(-0.45, 3.45)
        count_text = "\n".join((
            f"Elastic: {counts[4]}",
            f"Ionization: {counts[2]}",
            f"Ly-alpha: {counts[3]}",
            f"State change: {counts[1]}",
        ))
        axes[3].text(
            0.98, 0.97,
            f"Counts in {args.altitude_min:g}-{args.altitude_max:g} km\n"
            + count_text,
            transform=axes[3].transAxes,
            ha="right", va="top", fontsize=7,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "0.75",
                "alpha": 0.90,
            },
        )
    axes[4].set_ylabel("Scattering angle (deg)")
    axes[4].set_ylim(-3, 183)
    axes[4].set_yticks((0, 30, 60, 90, 120, 150, 180))
    axes[4].set_xlabel("Time (s)")
    axes[0].legend(ncol=2, fontsize=7, loc="best")

    # Lowercase bold panel labels follow common Nature-family conventions.
    for label, axis in zip(("a", "b", "c", "d", "e"), axes):
        axis.text(
            0.02, 0.96, label, transform=axis.transAxes,
            ha="left", va="top", fontweight="bold", fontsize=9,
        )
    for ax in axes:
        ax.grid(True, color="0.90", lw=0.55)
        ax.tick_params(direction="out", width=0.8, length=3)
    fig.suptitle(
        f"Single {species} trajectory, 400 km/s, atmospheric segment"
    )

    output = args.output or args.mat_file.with_suffix(".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    # Save one high-resolution PNG for the example and README preview.
    output_stem = output.with_suffix("")
    fig.savefig(output_stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"png={output_stem.with_suffix('.png').resolve()}")


if __name__ == "__main__":
    main()
