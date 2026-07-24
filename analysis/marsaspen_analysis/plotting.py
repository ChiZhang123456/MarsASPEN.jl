from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from .io import load_history_mat

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"

EVENT_NAMES = ("initial", "transport", "collision", "final")
TARGET_NAMES = ("none", "CO2", "O", "N2")
REACTION_NAMES = ("none", "state change", "ionization", "Ly-alpha", "elastic")
REACTION_COLORS = {
    1: "#d62728",
    2: "#ff7f0e",
    3: "#2ca02c",
    4: "#4d4d4d",
}


def one_dim(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    return np.asarray(data[key]).squeeze()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot detailed Aspen.jl MAT output.")
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    data = load_history_mat(args.mat_file)
    out = args.output_dir or args.mat_file.with_suffix("")
    out.mkdir(parents=True, exist_ok=True)

    pid = one_dim(data, "particle_id").astype(int)
    event_type = one_dim(data, "event_type").astype(int)
    time_s = one_dim(data, "time_s")
    altitude = one_dim(data, "altitude_km")
    energy = one_dim(data, "energy_ev")
    energy_before = one_dim(data, "energy_before_ev")
    charge = one_dim(data, "charge_state").astype(int)
    target = one_dim(data, "target_code").astype(int)
    reaction = one_dim(data, "reaction_code").astype(int)
    energy_loss = one_dim(data, "energy_loss_ev")
    x = one_dim(data, "x_m") / 1000.0
    y = one_dim(data, "y_m") / 1000.0
    z = one_dim(data, "z_m") / 1000.0
    vx = one_dim(data, "vx_m_s")
    vy = one_dim(data, "vy_m_s")
    vz = one_dim(data, "vz_m_s")
    vx_before = one_dim(data, "vx_before_m_s")
    vy_before = one_dim(data, "vy_before_m_s")
    vz_before = one_dim(data, "vz_before_m_s")
    speed = np.sqrt(vx**2 + vy**2 + vz**2) / 1000.0

    particle_ids = np.unique(pid)
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(particle_ids), 2)))
    fig, axes = plt.subplots(5, 1, figsize=(11, 14), sharex=True, constrained_layout=True)
    for color, particle_id in zip(colors, particle_ids):
        mask = pid == particle_id
        label = f"Particle {particle_id}"
        axes[0].plot(time_s[mask], altitude[mask], color=color, lw=0.9, label=label)
        axes[1].plot(time_s[mask], energy[mask], color=color, lw=0.9)
        axes[2].plot(time_s[mask], speed[mask], color=color, lw=0.9)
        axes[3].step(time_s[mask], charge[mask], where="post", color=color, lw=0.9)
        axes[4].plot(time_s[mask], x[mask], color=color, lw=0.9)

    collision = event_type == 2
    for code, color in REACTION_COLORS.items():
        mask = collision & (reaction == code)
        axes[0].scatter(time_s[mask], altitude[mask], s=11, color=color,
                        label=REACTION_NAMES[code], zorder=4)
        axes[1].scatter(time_s[mask], energy[mask], s=11, color=color, zorder=4)
    axes[0].set_ylabel("Altitude (km)")
    axes[1].set_ylabel("Energy (eV)")
    axes[2].set_ylabel("Speed (km/s)")
    axes[3].set_ylabel("Charge state")
    axes[3].set_yticks([0, 1], ["H", "H+"])
    axes[4].set_ylabel("MSO X (km)")
    axes[4].set_xlabel("Time (s)")
    axes[0].legend(ncol=3, fontsize=8)
    fig.suptitle(f"ASPEN Julia detailed trajectories, {len(particle_ids)} particles")
    overview = out / "trajectory_energy_reactions.png"
    fig.savefig(overview, dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    for color, particle_id in zip(colors, particle_ids):
        mask = pid == particle_id
        transverse = np.sqrt(y[mask] ** 2 + z[mask] ** 2)
        ax.plot(x[mask], transverse, color=color, lw=1.0, label=f"Particle {particle_id}")
    ax.set_xlabel("MSO X (km)")
    ax.set_ylabel(r"$\sqrt{Y^2+Z^2}$ (km)")
    ax.set_title("Particle trajectories")
    ax.legend(fontsize=8, ncol=2)
    trajectory = out / "trajectory_x_transverse.png"
    fig.savefig(trajectory, dpi=220)
    plt.close(fig)

    individual_dir = out / "individual_particles"
    individual_dir.mkdir(parents=True, exist_ok=True)
    individual_files = []
    for particle_id in particle_ids:
        particle = pid == particle_id
        particle_collision = particle & collision
        fig, ax = plt.subplots(figsize=(10, 6.4), constrained_layout=True)
        ax.plot(
            time_s[particle],
            altitude[particle],
            color="black",
            lw=1.5,
            zorder=1,
            label="trajectory",
        )
        for code in (4, 3, 2, 1):
            mask = particle_collision & (reaction == code)
            if not np.any(mask):
                continue
            marker_size = 18 if code == 4 else 38
            marker_alpha = 0.65 if code == 4 else 0.95
            ax.scatter(
                time_s[mask],
                altitude[mask],
                s=marker_size,
                color=REACTION_COLORS[code],
                edgecolors="none",
                alpha=marker_alpha,
                label=REACTION_NAMES[code],
                zorder=2 if code == 4 else 3,
            )
        ax.axhline(80.0, color="#e41a1c", ls="--", lw=1.2, label="lower boundary")
        ax.axhline(1000.0, color="#e41a1c", ls="--", lw=1.2, label="upper boundary")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Altitude (km)")
        ax.set_title(
            f"Particle {particle_id}, H-ENA at 600 km, "
            r"$\mathbf{V}_0=[-400,0,0]$ km/s"
        )
        ax.grid(True, color="0.85", lw=0.7, alpha=0.7)
        ax.legend(loc="best", fontsize=9, frameon=True)
        ax.set_ylim(60.0, 1040.0)
        filename = individual_dir / f"particle_{particle_id:02d}_altitude_reactions.png"
        fig.savefig(filename, dpi=220)
        plt.close(fig)
        individual_files.append(filename)

    collision_indices = np.flatnonzero(collision)
    csv_file = out / "chemical_reaction_events.csv"
    fields = [
        "particle_id", "time_s", "altitude_km", "target", "reaction",
        "energy_before_ev", "energy_after_ev", "energy_loss_ev",
        "speed_before_km_s", "speed_after_km_s",
        "vx_before_km_s", "vy_before_km_s", "vz_before_km_s",
        "vx_after_km_s", "vy_after_km_s", "vz_after_km_s", "charge_state_after",
    ]
    with csv_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i in collision_indices:
            writer.writerow({
                "particle_id": pid[i],
                "time_s": f"{time_s[i]:.9g}",
                "altitude_km": f"{altitude[i]:.9g}",
                "target": TARGET_NAMES[target[i]],
                "reaction": REACTION_NAMES[reaction[i]],
                "energy_before_ev": f"{energy_before[i]:.9g}",
                "energy_after_ev": f"{energy[i]:.9g}",
                "energy_loss_ev": f"{energy_loss[i]:.9g}",
                "speed_before_km_s": f"{np.sqrt(vx_before[i]**2 + vy_before[i]**2 + vz_before[i]**2) / 1000:.9g}",
                "speed_after_km_s": f"{speed[i]:.9g}",
                "vx_before_km_s": f"{vx_before[i] / 1000:.9g}",
                "vy_before_km_s": f"{vy_before[i] / 1000:.9g}",
                "vz_before_km_s": f"{vz_before[i] / 1000:.9g}",
                "vx_after_km_s": f"{vx[i] / 1000:.9g}",
                "vy_after_km_s": f"{vy[i] / 1000:.9g}",
                "vz_after_km_s": f"{vz[i] / 1000:.9g}",
                "charge_state_after": charge[i],
            })

    print(f"n_particles={len(particle_ids)}")
    print(f"n_collision_events={collision_indices.size}")
    for code in range(1, 5):
        print(f"n_{REACTION_NAMES[code].replace(' ', '_')}={np.count_nonzero(collision & (reaction == code))}")
    print(f"overview_png={overview.resolve()}")
    print(f"trajectory_png={trajectory.resolve()}")
    print(f"individual_particle_dir={individual_dir.resolve()}")
    print(f"n_individual_particle_plots={len(individual_files)}")
    print(f"reaction_csv={csv_file.resolve()}")


if __name__ == "__main__":
    main()
