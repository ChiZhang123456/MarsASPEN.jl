"""Plot local radial number-flux profiles for H ENA and H+.

The input MAT file is produced by either radial-flux Julia example. At each
spherical altitude crossing, Julia accumulates

    F_i = Wn_i * abs(Vr_i)

into a downward or upward array. ``Wn`` is in m^-3 and local radial velocity
``Vr = v dot r_hat`` is in m s^-1, so every profile is in m^-2 s^-1.

Panels a and b show positive directional magnitudes. Panel c shows the signed
outward radial flux

    F_r = F_upward - F_downward.

Thus negative values in panel c indicate net precipitation, while positive
values indicate net escape. Repeated upward and downward crossings contribute
to the directional panels, but cancel with the appropriate sign in the net
profile when their crossing speeds and charge states are unchanged.
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

from marsaspen_analysis import (  # noqa: E402
    load_history_mat,
    mat_string,
    radial_flux_from_mat,
)

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

SPECIES = (("H ENA", "#4C72B0"), (r"H$^+$", "#C44E52"))


def vector(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    """Return one MAT row or column as a one-dimensional float array."""
    return np.asarray(data[key], dtype=float).squeeze()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot local radial flux profiles from a MarsASPEN MAT file."
    )
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    data = load_history_mat(args.mat_file)
    profiles = radial_flux_from_mat(data)
    altitude = profiles["altitude_km"]
    downward = profiles["downward"]
    upward = profiles["upward"]
    signed_outward = profiles["signed_outward"]
    nominal_flux = float(vector(data, "nominal_bulk_flux_m2_s"))
    initial_species = mat_string(data["initial_species"]).replace("_", " ")
    source_label = r"H$^+$" if "plus" in initial_species.lower() else "H ENA"

    fig, axes = plt.subplots(
        1, 3, figsize=(7.2, 4.4), sharey=True, constrained_layout=True
    )
    directional_values = (downward, upward)
    directional_titles = ("Downward radial flux", "Upward radial flux")

    for axis, values, title in zip(
        axes[:2], directional_values, directional_titles
    ):
        for charge_index, (label, color) in enumerate(SPECIES):
            # Logarithmic axes cannot display zero. Masking preserves genuinely
            # empty altitude ranges instead of inventing a numerical floor.
            profile = np.ma.masked_less_equal(values[:, charge_index], 0.0)
            axis.plot(profile, altitude, color=color, lw=1.25, label=label)
        axis.set_xscale("log")
        axis.set_xlim(1.0e6, 3.0e12)
        axis.set_xlabel(r"Flux (m$^{-2}$ s$^{-1}$)")
        axis.set_title(title)
        axis.grid(True, which="major", color="0.90", lw=0.55)

    # The signed panel uses a linear scale. Negative values are net downward
    # precipitation and positive values are net outward escape.
    for charge_index, (label, color) in enumerate(SPECIES):
        axes[2].plot(
            signed_outward[:, charge_index], altitude,
            color=color, lw=1.25, label=label,
        )
    axes[2].axvline(0.0, color="0.35", lw=0.8)
    axes[2].set_xlim(-2.2e12, 2.2e12)
    axes[2].ticklabel_format(
        axis="x", style="sci", scilimits=(0, 0), useMathText=True
    )
    axes[2].set_xlabel(r"Signed $F_r$ (m$^{-2}$ s$^{-1}$)")
    axes[2].set_title("Net outward radial flux")
    axes[2].grid(True, which="major", color="0.90", lw=0.55)

    axes[0].set_ylabel("Altitude (km)")
    axes[0].set_ylim(80, 600)
    axes[0].legend(loc="best", fontsize=7)
    for label, axis in zip(("a", "b", "c"), axes):
        axis.text(
            0.03, 0.97, label, transform=axis.transAxes,
            ha="left", va="top", fontweight="bold", fontsize=9,
        )
    fig.suptitle(
        f"Local radial flux from 100,000 {source_label} source particles\n"
        r"600 km, 400 km/s, $T=10$ eV, $n=5$ cm$^{-3}$; "
        rf"$nU={nominal_flux:.2e}$ m$^{{-2}}$ s$^{{-1}}$",
        fontsize=9,
    )

    output = args.output or args.mat_file.with_name("radial_flux_profile.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={output.resolve()}")
    print(f"top_Hplus_downward_flux_m2_s={downward[-1, 1]:.9g}")
    print(f"top_HENA_upward_flux_m2_s={upward[-1, 0]:.9g}")
    print(f"top_Hplus_signed_outward_flux_m2_s={signed_outward[-1, 1]:.9g}")
    print(f"top_HENA_signed_outward_flux_m2_s={signed_outward[-1, 0]:.9g}")


if __name__ == "__main__":
    main()
