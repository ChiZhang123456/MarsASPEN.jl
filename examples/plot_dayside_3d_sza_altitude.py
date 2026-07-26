"""Plot SZA-altitude diagnostics from the full 3D Monte Carlo output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm

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

    fields = [
        np.asarray(moments["total_number_density_m3"], dtype=float),
        np.asarray(moments["total_flux_m2_s"], dtype=float),
        np.asarray(moments["downward_radial_flux_m2_s"], dtype=float),
        np.asarray(moments["upward_radial_flux_m2_s"], dtype=float),
        np.asarray(moments["signed_radial_flux_m2_s"], dtype=float),
        np.asarray(reactions["reaction_rate_by_channel_m3_s1"], dtype=float)[1],
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
        "H + H$^+$ number density",
        "Total scalar flux",
        "Downward radial flux",
        "Upward radial flux",
        "Signed outward radial flux",
        "Total target ionization rate",
        "H Ly-alpha volume emission rate",
        "Projectile energy transfer rate",
    )
    colorbar_labels = (
        r"Density (m$^{-3}$)",
        r"Flux (m$^{-2}$ s$^{-1}$)",
        r"Downward flux (m$^{-2}$ s$^{-1}$)",
        r"Upward flux (m$^{-2}$ s$^{-1}$)",
        r"Signed flux (m$^{-2}$ s$^{-1}$)",
        r"Ionization rate (m$^{-3}$ s$^{-1}$)",
        r"VER (photons m$^{-3}$ s$^{-1}$)",
        r"Energy transfer (W m$^{-3}$)",
    )

    fig, axes = plt.subplots(
        2, 4, figsize=(10.0, 5.3), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for panel, (axis, values, title, colorbar_label) in enumerate(
        zip(axes.flat, sza_altitude, titles, colorbar_labels)
    ):
        if panel == 4:
            maximum = np.nanmax(np.abs(values))
            norm = SymLogNorm(
                linthresh=max(maximum / 1000.0, 1.0),
                vmin=-maximum, vmax=maximum, base=10,
            )
            plotted = values
            cmap = "coolwarm"
        else:
            low, high = positive_log_limits(values)
            norm = LogNorm(vmin=low, vmax=high)
            plotted = np.ma.masked_less_equal(values, 0.0)
            cmap = "turbo"
        mesh = axis.pcolormesh(
            SZA_EDGES_DEG, displayed_altitude_edges, plotted,
            shading="flat", cmap=cmap, norm=norm, rasterized=True,
        )
        axis.set(
            xlim=(0, 180), ylim=(ALTITUDE_MIN_KM, ALTITUDE_MAX_KM),
            title=title,
        )
        if panel >= 4:
            axis.set_xlabel("Solar zenith angle (deg)")
        if panel % 4 == 0:
            axis.set_ylabel("Altitude (km)")
        axis.axvline(90, color="0.25", lw=0.7, ls=":", alpha=0.9)
        axis.text(
            0.025, 0.975, "abcdefgh"[panel],
            transform=axis.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
        )
        colorbar = fig.colorbar(mesh, ax=axis, pad=0.015, fraction=0.048)
        colorbar.set_label(colorbar_label)

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
