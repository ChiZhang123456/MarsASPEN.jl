"""Plot SZA-altitude diagnostics from the full 3D Monte Carlo output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis import load_history_mat  # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

ALTITUDE_MIN_KM = 100.0
ALTITUDE_MAX_KM = 300.0
SZA_EDGES_DEG = np.arange(0.0, 180.0 + 5.0, 5.0)


def sza_area_mean(
    values: np.ndarray,
    latitude_edges_deg: np.ndarray,
    longitude_centers_deg: np.ndarray,
) -> np.ndarray:
    """Area-average a lon-lat-alt field within each 5 degree SZA annulus."""
    latitude_centers = 0.5 * (
        latitude_edges_deg[:-1] + latitude_edges_deg[1:]
    )
    latitude_weight = np.diff(
        np.sin(np.deg2rad(latitude_edges_deg))
    )
    longitude_rad = np.deg2rad(longitude_centers_deg)
    latitude_rad = np.deg2rad(latitude_centers)
    sza = np.rad2deg(np.arccos(np.clip(
        np.cos(latitude_rad)[:, None] *
        np.cos(longitude_rad)[None, :],
        -1.0, 1.0,
    )))
    cell_weight = np.broadcast_to(
        latitude_weight[:, None], sza.shape
    )
    output = np.zeros((values.shape[0], SZA_EDGES_DEG.size - 1))
    for index in range(SZA_EDGES_DEG.size - 1):
        selected = (
            (sza >= SZA_EDGES_DEG[index]) &
            (sza < SZA_EDGES_DEG[index + 1])
        )
        weight = cell_weight[selected]
        if weight.size:
            output[:, index] = (
                np.sum(values[:, selected] * weight[None, :], axis=1) /
                np.sum(weight)
            )
    return output


def positive_log_limits(values: np.ndarray) -> tuple[float, float]:
    """Choose robust logarithmic limits for one SZA-altitude panel."""
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        return 1.0, 10.0
    low = 10.0 ** np.floor(np.log10(np.percentile(positive, 3)))
    high = 10.0 ** np.ceil(np.log10(np.max(positive)))
    return low, max(high, low * 10.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "prefix", nargs="?", type=Path,
        default=REPO / "examples" / "output" /
        "dayside_hplus_10000000_3d",
    )
    parser.add_argument(
        "--output", type=Path,
        default=REPO / "examples" / "figures" /
        "dayside_hplus_10000000_3d_sza_altitude_8panel.png",
    )
    args = parser.parse_args()
    prefix = str(args.prefix)
    moments = load_history_mat(Path(prefix + "_moments.mat"))
    reactions = load_history_mat(Path(prefix + "_reactions.mat"))
    energy = load_history_mat(Path(prefix + "_energy.mat"))

    altitude_centers = np.ravel(moments["altitude_centers_km"])
    altitude_edges = np.ravel(moments["altitude_edges_km"])
    latitude_edges = np.ravel(moments["latitude_edges_deg"])
    longitude_centers = np.ravel(moments["longitude_centers_deg"])
    altitude_selected = (
        (altitude_centers >= ALTITUDE_MIN_KM) &
        (altitude_centers < ALTITUDE_MAX_KM)
    )
    first = int(np.flatnonzero(altitude_selected)[0])
    last = int(np.flatnonzero(altitude_selected)[-1])
    displayed_altitude_edges = altitude_edges[first:last + 2]

    number_density = np.asarray(
        moments["number_density_by_charge_m3"], dtype=float
    )
    downward_flux = np.asarray(
        moments["downward_radial_flux_by_charge_m2_s"], dtype=float
    )
    upward_flux = np.asarray(
        moments["upward_radial_flux_by_charge_m2_s"], dtype=float
    )
    ionization = np.asarray(
        reactions["ionization_rate_by_target_m3_s1"], dtype=float
    )
    fields = [
        number_density[1],
        number_density[0],
        downward_flux[1],
        upward_flux[1],
        downward_flux[0],
        upward_flux[0],
        ionization[1],
        ionization[0],
        np.asarray(
            reactions["total_lya_volume_emission_rate_photons_m3_s1"],
            dtype=float,
        ),
        np.asarray(energy["total_energy_transfer_w_m3"], dtype=float),
    ]
    sza_altitude = [
        sza_area_mean(values, latitude_edges, longitude_centers)[
            altitude_selected
        ]
        for values in fields
    ]
    titles = (
        "H$^+$ number density",
        "H-ENA number density",
        "H$^+$ downward radial flux",
        "H$^+$ upward radial flux",
        "H-ENA downward radial flux",
        "H-ENA upward radial flux",
        "O ionization rate",
        "CO$_2$ ionization rate",
        "H Ly-alpha volume emission rate",
        "Energy deposition rate",
    )
    row_colorbar_labels = (
        r"Density (m$^{-3}$)",
        r"Flux (m$^{-2}$ s$^{-1}$)",
        r"Flux (m$^{-2}$ s$^{-1}$)",
        r"Ionization rate (m$^{-3}$ s$^{-1}$)",
    )

    fig, axes = plt.subplots(
        5, 2, figsize=(7.2, 10.0), sharex=True, sharey=True,
        layout="constrained",
    )
    meshes = []
    for panel, (axis, values, title) in enumerate(
        zip(axes.flat, sza_altitude, titles)
    ):
        row = panel // 2
        if panel % 2 == 0 and row < 4:
            paired = np.concatenate((
                sza_altitude[panel].ravel(),
                sza_altitude[panel + 1].ravel(),
            ))
            low, high = positive_log_limits(paired)
            shared_norm = LogNorm(vmin=low, vmax=high)
        if row < 4:
            norm = shared_norm
        else:
            low, high = positive_log_limits(values)
            norm = LogNorm(vmin=low, vmax=high)
        plotted = np.ma.masked_less_equal(values, 0.0)
        mesh = axis.pcolormesh(
            SZA_EDGES_DEG, displayed_altitude_edges, plotted,
            shading="flat", cmap="turbo", norm=norm, rasterized=True,
        )
        meshes.append(mesh)
        axis.set(
            xlim=(0, 180), ylim=(ALTITUDE_MIN_KM, ALTITUDE_MAX_KM),
            title=title,
        )
        if row == 4:
            axis.set_xlabel("Solar zenith angle (deg)")
        if panel % 2 == 0:
            axis.set_ylabel("Altitude (km)")
        axis.axvline(90, color="0.25", lw=0.7, ls=":", alpha=0.9)
        axis.text(
            0.025, 0.975, "abcdefghij"[panel],
            transform=axis.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
        )

    for row, colorbar_label in enumerate(row_colorbar_labels):
        colorbar = fig.colorbar(
            meshes[2 * row], ax=axes[row, :],
            pad=0.02, fraction=0.045, aspect=32,
        )
        colorbar.set_label(colorbar_label)
    colorbar = fig.colorbar(
        meshes[8], ax=axes[4, 0], pad=0.02, fraction=0.045, aspect=32,
    )
    colorbar.set_label(r"VER (photons m$^{-3}$ s$^{-1}$)")
    colorbar = fig.colorbar(
        meshes[9], ax=axes[4, 1], pad=0.02, fraction=0.045, aspect=32,
    )
    colorbar.set_label(r"Energy deposition rate (W m$^{-3}$)")

    fig.suptitle(
        "Uniform-dayside H$^+$ Monte Carlo SZA-altitude diagnostics\n"
        r"10,000,000 particles, $\mathbf{U}=(-400,0,0)$ km/s, "
        r"$T=10$ eV, $n=5$ cm$^{-3}$, "
        r"5$^\circ$ SZA $\times$ 1 km altitude bins",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={args.output.resolve()}")
    for title, values in zip(titles, sza_altitude):
        print(f"{title}: max_abs={np.nanmax(np.abs(values)):.9g}")


if __name__ == "__main__":
    main()
