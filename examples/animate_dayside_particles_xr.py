"""Animate MarsASPEN particle positions in the MSO X-R plane.

R is the cylindrical distance from the MSO X axis:

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
                "r": np.hypot(y_km[indices], z_km[indices]),
                "charge": charge[indices],
            }
        )
    return tracks


def _sample_tracks(
    tracks: list[dict[str, np.ndarray]], frame_times: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate positions and sample charge states at every frame."""
    n_frames = frame_times.size
    n_particles = len(tracks)
    x = np.full((n_frames, n_particles), np.nan)
    r = np.full_like(x, np.nan)
    charge = np.full((n_frames, n_particles), -1, dtype=np.int8)
    active = np.zeros((n_frames, n_particles), dtype=bool)

    for column, track in enumerate(tracks):
        times = track["time"]
        mask = (frame_times >= times[0]) & (frame_times <= times[-1])
        active[mask, column] = True
        x[mask, column] = np.interp(frame_times[mask], times, track["x"])
        r[mask, column] = np.interp(frame_times[mask], times, track["r"])
        state_index = np.searchsorted(times, frame_times[mask], side="right") - 1
        charge[mask, column] = track["charge"][state_index]
    return x, r, charge, active


def make_animation(
    input_file: Path,
    output_file: Path,
    time_step_s: float = 0.1,
    playback_ms: int = 120,
) -> None:
    """Render the X-R animation as an animated GIF."""
    data = load_history_mat(input_file)
    tracks = _particle_tracks(data)
    maximum_time = max(track["time"][-1] for track in tracks)
    frame_times = np.arange(0.0, maximum_time + 0.5 * time_step_s, time_step_s)
    x, r, charge, active = _sample_tracks(tracks, frame_times)

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

    figure, axis = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)
    angle = np.linspace(0.0, np.pi, 500)
    mars_x = MARS_RADIUS_KM * np.cos(angle)
    mars_r = MARS_RADIUS_KM * np.sin(angle)
    injection_radius = MARS_RADIUS_KM + INJECTION_ALTITUDE_KM
    injection_x = injection_radius * np.cos(angle)
    injection_r = injection_radius * np.sin(angle)

    axis.fill_between(
        mars_x, 0.0, mars_r, color="#C86B4A", alpha=0.85, linewidth=0
    )
    axis.plot(mars_x, mars_r, color="#7A3827", linewidth=0.9)
    axis.plot(
        injection_x,
        injection_r,
        color="#7A7A7A",
        linewidth=0.8,
        linestyle=(0, (3, 2)),
    )
    hplus_points = axis.scatter(
        [], [], s=8, color="#2878B5", alpha=0.72, linewidths=0, zorder=4
    )
    hena_points = axis.scatter(
        [], [], s=8, color="#E07B39", alpha=0.78, linewidths=0, zorder=5
    )
    time_text = axis.text(
        0.025,
        0.96,
        "",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10,
    )
    count_text = axis.text(
        0.025,
        0.895,
        "",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="#4D4D4D",
    )

    axis.set_xlim(-4_450, 4_450)
    axis.set_ylim(0, 4_450)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("MSO X (km)")
    axis.set_ylabel(r"$R=\sqrt{Y^2+Z^2}$ (km)")
    axis.set_title("Dayside solar-wind proton transport")
    axis.legend(
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
        loc="upper right",
        frameon=False,
    )

    empty = np.empty((0, 2))

    def update(frame: int):
        valid = active[frame]
        proton = valid & (charge[frame] == 1)
        neutral = valid & (charge[frame] == 0)
        hplus_points.set_offsets(
            np.column_stack((x[frame, proton], r[frame, proton]))
            if np.any(proton)
            else empty
        )
        hena_points.set_offsets(
            np.column_stack((x[frame, neutral], r[frame, neutral]))
            if np.any(neutral)
            else empty
        )
        time_text.set_text(f"Simulation time = {frame_times[frame]:.1f} s")
        count_text.set_text(
            f"Active particles: {np.count_nonzero(valid):,}   "
            f"H+: {np.count_nonzero(proton):,}   "
            f"H ENA: {np.count_nonzero(neutral):,}"
        )
        return hplus_points, hena_points, time_text, count_text

    animation = FuncAnimation(
        figure,
        update,
        frames=frame_times.size,
        interval=playback_ms,
        blit=True,
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    animation.save(
        output_file,
        writer=PillowWriter(fps=max(1, round(1000 / playback_ms))),
        dpi=150,
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
