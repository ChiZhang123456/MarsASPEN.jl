"""Animate MarsASPEN particle positions in four MSO projections.

The panels are X-R, X-Y, X-Z, and Y-Z, where R is the cylindrical distance
from the MSO X axis:

    R = sqrt(Y**2 + Z**2)

The detailed, irregularly sampled particle histories are linearly
interpolated onto a common 0.1 s time grid. Charge state is held constant
between recorded events, so charge exchange appears at the correct event
time.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis.io import load_history_mat


MARS_RADIUS_KM = 3389.5
INJECTION_ALTITUDE_KM = 600.0


def _vector(data: dict[str, np.ndarray], key: str, dtype=float) -> np.ndarray:
    """Return one MAT variable as a flat NumPy vector."""
    return np.asarray(data[key], dtype=dtype).reshape(-1)


def _particle_tracks(data: dict[str, np.ndarray]) -> list[dict[str, np.ndarray]]:
    """Split flattened event columns into monotonic per-particle tracks."""
    particle_id = _vector(data, "particle_id", int)
    time_s = _vector(data, "time_s")
    x_km = _vector(data, "x_m") / 1000.0
    y_km = _vector(data, "y_m") / 1000.0
    z_km = _vector(data, "z_m") / 1000.0
    charge = _vector(data, "charge_state", int)

    tracks: list[dict[str, np.ndarray]] = []
    for pid in np.unique(particle_id):
        indices = np.flatnonzero(particle_id == pid)
        order = np.argsort(time_s[indices], kind="stable")
        indices = indices[order]
        times = time_s[indices]

        # A collision row and its preceding transport row can have the same
        # time and position. Retain the last row so the updated charge state
        # is used immediately after the collision.
        keep = np.r_[np.flatnonzero(np.diff(times) > 0), times.size - 1]
        indices = indices[keep]
        tracks.append(
            {
                "time": time_s[indices],
                "x": x_km[indices],
                "y": y_km[indices],
                "z": z_km[indices],
                "charge": charge[indices],
            }
        )
    return tracks


def _sample_tracks(
    tracks: list[dict[str, np.ndarray]], frame_times: np.ndarray
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Interpolate positions and sample charge states at every frame."""
    n_frames = frame_times.size
    n_particles = len(tracks)
    x = np.full((n_frames, n_particles), np.nan)
    y = np.full_like(x, np.nan)
    z = np.full_like(x, np.nan)
    charge = np.full((n_frames, n_particles), -1, dtype=np.int8)
    active = np.zeros((n_frames, n_particles), dtype=bool)

    for column, track in enumerate(tracks):
        times = track["time"]
        mask = (frame_times >= times[0]) & (frame_times <= times[-1])
        active[mask, column] = True
        x[mask, column] = np.interp(frame_times[mask], times, track["x"])
        y[mask, column] = np.interp(frame_times[mask], times, track["y"])
        z[mask, column] = np.interp(frame_times[mask], times, track["z"])
        state_index = np.searchsorted(times, frame_times[mask], side="right") - 1
        charge[mask, column] = track["charge"][state_index]
    r = np.hypot(y, z)
    return x, y, z, r, charge, active


def make_animation(
    input_file: Path,
    output_file: Path,
    time_step_s: float = 0.1,
    playback_ms: int = 120,
) -> None:
    """Render four synchronized particle projections as an animated GIF."""
    data = load_history_mat(input_file)
    tracks = _particle_tracks(data)
    maximum_time = max(track["time"][-1] for track in tracks)
    frame_times = np.arange(0.0, maximum_time + 0.5 * time_step_s, time_step_s)
    x, y, z, r, charge, active = _sample_tracks(tracks, frame_times)

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )

    figure, axes = plt.subplots(
        2, 2, figsize=(8.6, 8.0), constrained_layout=True
    )
    angle_half = np.linspace(0.0, np.pi, 500)
    angle_full = np.linspace(0.0, 2.0 * np.pi, 700)
    mars_x = MARS_RADIUS_KM * np.cos(angle_half)
    mars_r = MARS_RADIUS_KM * np.sin(angle_half)
    injection_radius = MARS_RADIUS_KM + INJECTION_ALTITUDE_KM
    injection_x = injection_radius * np.cos(angle_half)
    injection_r = injection_radius * np.sin(angle_half)

    axes[0, 0].fill_between(
        mars_x, 0.0, mars_r, color="#C86B4A", alpha=0.85, linewidth=0
    )
    axes[0, 0].plot(mars_x, mars_r, color="#7A3827", linewidth=0.9)
    axes[0, 0].plot(
        injection_x,
        injection_r,
        color="#7A7A7A",
        linewidth=0.8,
        linestyle=(0, (3, 2)),
    )

    # The three Cartesian projections show Mars as a full circular disk and
    # the 600 km source radius as a dashed circle.
    circle_x = MARS_RADIUS_KM * np.cos(angle_full)
    circle_y = MARS_RADIUS_KM * np.sin(angle_full)
    source_x = injection_radius * np.cos(angle_full)
    source_y = injection_radius * np.sin(angle_full)
    for axis in (axes[0, 1], axes[1, 0], axes[1, 1]):
        axis.fill(circle_x, circle_y, color="#C86B4A", alpha=0.85, linewidth=0)
        axis.plot(circle_x, circle_y, color="#7A3827", linewidth=0.9)
        axis.plot(
            source_x,
            source_y,
            color="#7A7A7A",
            linewidth=0.8,
            linestyle=(0, (3, 2)),
        )

    panel_coordinates = (
        (x, r),
        (x, y),
        (x, z),
        (y, z),
    )
    hplus_points = []
    hena_points = []
    for axis in axes.flat:
        hplus_points.append(
            axis.scatter(
                [], [], s=6, color="#2878B5", alpha=0.70,
                linewidths=0, zorder=4
            )
        )
        hena_points.append(
            axis.scatter(
                [], [], s=6, color="#E07B39", alpha=0.78,
                linewidths=0, zorder=5
            )
        )

    axis_limit = 4_450
    for axis in axes.flat:
        axis.set_xlim(-axis_limit, axis_limit)
        axis.set_ylim(-axis_limit, axis_limit)
        axis.set_aspect("equal", adjustable="box")
    axes[0, 0].set_ylim(0, axis_limit)
    axes[0, 0].set_xlabel("MSO X (km)")
    axes[0, 0].set_ylabel(r"$R=\sqrt{Y^2+Z^2}$ (km)")
    axes[0, 0].set_title("X-R")
    axes[0, 1].set_xlabel("MSO X (km)")
    axes[0, 1].set_ylabel("MSO Y (km)")
    axes[0, 1].set_title("X-Y")
    axes[1, 0].set_xlabel("MSO X (km)")
    axes[1, 0].set_ylabel("MSO Z (km)")
    axes[1, 0].set_title("X-Z")
    axes[1, 1].set_xlabel("MSO Y (km)")
    axes[1, 1].set_ylabel("MSO Z (km)")
    axes[1, 1].set_title("Y-Z")

    title_text = figure.suptitle(
        "", fontsize=11, linespacing=1.45, color="#222222"
    )
    axes[0, 0].legend(
        handles=[
            Line2D(
                [0], [0], marker="o", linestyle="none", color="#2878B5",
                markersize=5, label="H+"
            ),
            Line2D(
                [0], [0], marker="o", linestyle="none", color="#E07B39",
                markersize=5, label="H ENA"
            ),
            Line2D(
                [0], [0], color="#7A7A7A", linestyle=(0, (3, 2)),
                linewidth=0.8, label="600 km injection surface"
            ),
        ],
        loc="upper left",
        fontsize=8,
        frameon=False,
    )

    empty = np.empty((0, 2))

    def update(frame: int):
        valid = active[frame]
        proton = valid & (charge[frame] == 1)
        neutral = valid & (charge[frame] == 0)
        for panel, (horizontal, vertical) in enumerate(panel_coordinates):
            hplus_points[panel].set_offsets(
                np.column_stack(
                    (horizontal[frame, proton], vertical[frame, proton])
                )
                if np.any(proton)
                else empty
            )
            hena_points[panel].set_offsets(
                np.column_stack(
                    (horizontal[frame, neutral], vertical[frame, neutral])
                )
                if np.any(neutral)
                else empty
            )
        title_text.set_text(
            "Dayside solar-wind proton transport\n"
            f"Time = {frame_times[frame]:.1f} s   "
            f"Active: {np.count_nonzero(valid):,}   "
            f"H+: {np.count_nonzero(proton):,}   "
            f"H ENA: {np.count_nonzero(neutral):,}"
        )
        return (*hplus_points, *hena_points, title_text)

    animation = FuncAnimation(
        figure,
        update,
        frames=frame_times.size,
        interval=playback_ms,
        blit=False,
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    animation.save(
        output_file,
        writer=PillowWriter(fps=max(1, round(1000 / playback_ms))),
        dpi=125,
    )
    plt.close(figure)
    print(f"frames={frame_times.size}")
    print(f"simulation_time_s={maximum_time:.6f}")
    print(f"output={output_file.resolve()}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=REPO / "examples" / "output" / "dayside_particle_animation_1000.mat",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "examples" / "figures" / "dayside_particles_xr.gif",
    )
    parser.add_argument("--time-step", type=float, default=0.1)
    parser.add_argument("--playback-ms", type=int, default=120)
    args = parser.parse_args()
    make_animation(
        args.input,
        args.output,
        time_step_s=args.time_step,
        playback_ms=args.playback_ms,
    )


if __name__ == "__main__":
    main()
