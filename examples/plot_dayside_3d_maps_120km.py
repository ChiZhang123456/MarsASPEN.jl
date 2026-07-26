"""Plot 120 to 121 km maps from the full 3D dayside Monte Carlo output."""

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


def positive_log_limits(values: np.ndarray) -> tuple[float, float]:
    """Return robust decade limits while retaining the footprint outskirts."""
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        return 1.0, 10.0
    low = 10.0 ** np.floor(np.log10(np.percentile(positive, 5)))
    high = 10.0 ** np.ceil(np.log10(np.max(positive)))
    return low, max(high, low * 10.0)


def longitude_shift(values: np.ndarray) -> np.ndarray:
    """Move the 0 degree MSO subsolar longitude to the map center."""
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
        np.asarray(moments["number_density_by_charge_m3"])[
            1, altitude_index
        ],
        np.asarray(moments["number_density_by_charge_m3"])[
            0, altitude_index
        ],
        np.asarray(moments["downward_radial_flux_by_charge_m2_s"])[
            1, altitude_index
        ],
        np.asarray(moments["upward_radial_flux_by_charge_m2_s"])[
            1, altitude_index
        ],
        np.asarray(reactions["ionization_rate_by_target_m3_s1"])[
            1, altitude_index
        ],
        np.asarray(reactions["ionization_rate_by_target_m3_s1"])[
            0, altitude_index
        ],
        np.asarray(
            reactions["total_lya_volume_emission_rate_photons_m3_s1"]
        )[altitude_index],
        np.asarray(energy["total_energy_transfer_w_m3"])[altitude_index],
    ]
    maps = [longitude_shift(values) for values in maps]
    titles = (
        "H$^+$ number density",
        "H-ENA number density",
        "H$^+$ downward radial flux",
        "H$^+$ upward radial flux",
        "O ionization rate",
        "CO$_2$ ionization rate",
        "H Ly-alpha volume emission rate",
        "Energy deposition rate",
    )

    fig, axes = plt.subplots(
        4, 2, figsize=(7.2, 8.2), sharex=True, sharey=True,
        constrained_layout=True,
    )
    shared_groups = ((0, 1), (2, 3), (4, 5))
    norms: dict[int, LogNorm] = {}
    for first, second in shared_groups:
        low, high = positive_log_limits(
            np.concatenate((maps[first].ravel(), maps[second].ravel()))
        )
        norm = LogNorm(vmin=low, vmax=high)
        norms[first] = norm
        norms[second] = norm
    for panel in (6, 7):
        low, high = positive_log_limits(maps[panel])
        norms[panel] = LogNorm(vmin=low, vmax=high)

    meshes = []
    for panel, (axis, values, title) in enumerate(
        zip(axes.flat, maps, titles)
    ):
        plotted = np.ma.masked_less_equal(values, 0.0)
        mesh = axis.pcolormesh(
            longitude_edges, latitude_edges, plotted,
            shading="flat", cmap="turbo", norm=norms[panel],
            rasterized=True,
        )
        meshes.append(mesh)
        axis.set(xlim=(-180, 180), ylim=(-90, 90), title=title)
        if panel >= 6:
            axis.set_xlabel("MSO longitude (deg)")
        if panel % 2 == 0:
            axis.set_ylabel("MSO latitude (deg)")
        axis.axvline(-90, color="0.25", lw=0.55, ls=":", alpha=0.8)
        axis.axvline(90, color="0.25", lw=0.55, ls=":", alpha=0.8)
        axis.set_xticks((-180, -90, 0, 90, 180))
        axis.text(
            0.02, 0.97, "abcdefgh"[panel],
            transform=axis.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
        )

    shared_colorbars = (
        ((0, 1), r"Number density (m$^{-3}$)"),
        ((2, 3), r"Radial flux (m$^{-2}$ s$^{-1}$)"),
        ((4, 5), r"Ionization rate (m$^{-3}$ s$^{-1}$)"),
    )
    for (first, second), label in shared_colorbars:
        colorbar = fig.colorbar(
            meshes[first],
            ax=(axes.flat[first], axes.flat[second]),
            pad=0.015, fraction=0.035,
        )
        colorbar.set_label(label)
    lya_colorbar = fig.colorbar(
        meshes[6], ax=axes.flat[6], pad=0.015, fraction=0.048
    )
    lya_colorbar.set_label(r"VER (photons m$^{-3}$ s$^{-1}$)")
    energy_colorbar = fig.colorbar(
        meshes[7], ax=axes.flat[7], pad=0.015, fraction=0.048
    )
    energy_colorbar.set_label(r"Energy deposition rate (W m$^{-3}$)")

    fig.suptitle(
        "Uniform-dayside H$^+$ Monte Carlo diagnostics at 120 km\n"
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
