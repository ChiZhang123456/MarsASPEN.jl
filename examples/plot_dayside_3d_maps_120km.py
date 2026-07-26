"""Plot 120 to 121 km maps from the full 3D dayside Monte Carlo output."""

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


def positive_log_limits(values: np.ndarray) -> tuple[float, float]:
    """Return robust decade limits while retaining the footprint outskirts."""
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        return 1.0, 10.0
    low = 10.0 ** np.floor(np.log10(np.percentile(positive, 5)))
    high = 10.0 ** np.ceil(np.log10(np.max(positive)))
    return low, max(high, low * 10.0)


def longitude_shift(values: np.ndarray) -> np.ndarray:
    """Convert the final array axis from 0..360 to -180..180 longitude."""
    return np.roll(values, -(values.shape[-1] // 2), axis=-1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "prefix", nargs="?", type=Path,
        default=REPO / "examples" / "output" /
        "dayside_hplus_10000000_3d",
    )
    parser.add_argument("--altitude", type=float, default=120.5)
    parser.add_argument(
        "--output", type=Path,
        default=REPO / "examples" / "figures" /
        "dayside_hplus_10000000_3d_maps_120km_8panel.png",
    )
    args = parser.parse_args()
    prefix = str(args.prefix)
    moments = load_history_mat(Path(prefix + "_moments.mat"))
    reactions = load_history_mat(Path(prefix + "_reactions.mat"))
    energy = load_history_mat(Path(prefix + "_energy.mat"))

    altitude_centers = np.ravel(moments["altitude_centers_km"])
    altitude_index = int(np.argmin(np.abs(altitude_centers - args.altitude)))
    altitude_center = altitude_centers[altitude_index]
    altitude_edges = np.ravel(moments["altitude_edges_km"])
    altitude_low = altitude_edges[altitude_index]
    altitude_high = altitude_edges[altitude_index + 1]
    latitude_edges = np.ravel(moments["latitude_edges_deg"])
    longitude_edges = np.linspace(-180.0, 180.0, 73)

    maps = [
        np.asarray(moments["total_number_density_m3"])[altitude_index],
        np.asarray(moments["total_flux_m2_s"])[altitude_index],
        np.asarray(moments["downward_radial_flux_m2_s"])[altitude_index],
        np.asarray(moments["upward_radial_flux_m2_s"])[altitude_index],
        np.asarray(moments["signed_radial_flux_m2_s"])[altitude_index],
        np.asarray(reactions["reaction_rate_by_channel_m3_s1"])[
            1, altitude_index
        ],
        np.asarray(
            reactions["total_lya_volume_emission_rate_photons_m3_s1"]
        )[altitude_index],
        np.asarray(energy["total_energy_transfer_w_m3"])[altitude_index],
    ]
    maps = [longitude_shift(values) for values in maps]
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
        zip(axes.flat, maps, titles, colorbar_labels)
    ):
        if panel == 4:
            maximum = np.nanmax(np.abs(values))
            norm = SymLogNorm(
                linthresh=max(maximum / 1000.0, 1.0),
                vmin=-maximum, vmax=maximum, base=10,
            )
            cmap = "coolwarm"
            plotted = values
        else:
            low, high = positive_log_limits(values)
            norm = LogNorm(vmin=low, vmax=high)
            cmap = "turbo"
            plotted = np.ma.masked_less_equal(values, 0.0)
        mesh = axis.pcolormesh(
            longitude_edges, latitude_edges, plotted,
            shading="flat", cmap=cmap, norm=norm, rasterized=True,
        )
        axis.set(xlim=(-180, 180), ylim=(-90, 90), title=title)
        if panel >= 4:
            axis.set_xlabel("MSO longitude (deg)")
        if panel % 4 == 0:
            axis.set_ylabel("MSO latitude (deg)")
        axis.axvline(-90, color="0.25", lw=0.55, ls=":", alpha=0.8)
        axis.axvline(90, color="0.25", lw=0.55, ls=":", alpha=0.8)
        axis.text(
            0.02, 0.97, "abcdefgh"[panel],
            transform=axis.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
        )
        colorbar = fig.colorbar(mesh, ax=axis, pad=0.015, fraction=0.048)
        if panel == 4:
            colorbar.set_ticks((-1.0e12, -1.0e9, 0.0, 1.0e9, 1.0e12))
        colorbar.set_label(colorbar_label)

    fig.suptitle(
        "Uniform-dayside H$^+$ Monte Carlo diagnostics\n"
        f"{altitude_low:.0f} to {altitude_high:.0f} km "
        f"(center {altitude_center:.1f} km), "
        r"10,000,000 particles, $\mathbf{U}=(-400,0,0)$ km/s, "
        r"$T=10$ eV, $n=5$ cm$^{-3}$",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"altitude_bin_km=({altitude_low}, {altitude_high})")
    print(f"output={args.output.resolve()}")
    for title, values in zip(titles, maps):
        print(f"{title}: max_abs={np.nanmax(np.abs(values)):.9g}")


if __name__ == "__main__":
    main()
