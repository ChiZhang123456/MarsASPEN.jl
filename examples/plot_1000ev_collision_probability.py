"""Calculate 1 km collision probabilities for a fixed 1000 eV projectile.

For every packaged season and solar-activity case, this script evaluates

    alpha(z) = sum_target n_target(z) sum_reaction sigma(E)

for H ENA and H+ at E = 1000 eV. The local probability over one downward
kilometer is 1 - exp[-alpha(z) * 1000 m]. The cumulative probability is
1 - exp[-integral_z^1000km alpha(s) ds].

This is a fixed-energy diagnostic. It does not reduce the projectile energy
after a collision and therefore does not replace a Monte Carlo trajectory.
It isolates how the atmosphere and cross sections set the first-collision
optical depth under different seasons and solar-activity levels.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

REPO = Path(__file__).resolve().parents[1]
ATMOSPHERE = REPO / "data" / "atmosphere"
CROSS_SECTIONS = REPO / "data" / "cross_sections"
FIGURES = REPO / "examples" / "figures"

LS_VALUES = (0, 90, 180, 270)
F107_VALUES = (70, 130, 200)
TARGETS = ("CO2", "O", "N2")
PROJECTILES = (("H", "H ENA", "#4C72B0"), ("Hplus", r"H$^+$", "#C44E52"))
ENERGY_EV = 1000.0
AMU_KG = 1.66053906660e-27
KB_J_K = 1.380649e-23
MARS_G0_M_S2 = 3.71
MARS_RADIUS_KM = 3388.25
TARGET_MASS_AMU = (44.01, 15.999, 28.014)

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def horizontal_profile(
    data: dict[str, np.ndarray], field: str, logarithmic: bool,
) -> np.ndarray:
    """Match Julia's periodic-longitude bilinear interpolation at 0, 0."""
    lon = np.asarray(data["lon_deg"], dtype=float).squeeze()
    lat = np.asarray(data["lat_deg"], dtype=float).squeeze()
    i0, i1 = lon.size - 1, 0
    wx = (360.0 - lon[i0]) / (lon[i1] + 360.0 - lon[i0])
    j0 = int(np.searchsorted(lat, 0.0) - 1)
    j1 = j0 + 1
    wy = (0.0 - lat[j0]) / (lat[j1] - lat[j0])
    values = np.asarray(data[field], dtype=float)
    if logarithmic:
        values = np.log(np.maximum(values, 1.0e-300))
    c0 = (1.0 - wx) * values[i0, j0, :] + wx * values[i1, j0, :]
    c1 = (1.0 - wx) * values[i0, j1, :] + wx * values[i1, j1, :]
    profile = (1.0 - wy) * c0 + wy * c1
    return np.exp(profile) if logarithmic else profile


def cold_model_profile(
    native_altitude: np.ndarray,
    native_density: np.ndarray,
    top_temperature_k: float,
    mass_amu: float,
    requested_altitude: np.ndarray,
) -> np.ndarray:
    """Match log interpolation, lower extrapolation, and hydrostatic top."""
    logn = np.log(np.maximum(native_density, 1.0e-300))
    result = np.exp(np.interp(
        np.minimum(requested_altitude, native_altitude[-1]),
        native_altitude, logn,
    ))
    below = requested_altitude < native_altitude[0]
    slope = (logn[1] - logn[0]) / (
        native_altitude[1] - native_altitude[0]
    )
    result[below] = np.exp(
        logn[0] + slope * (requested_altitude[below] - native_altitude[0])
    )
    above = requested_altitude > native_altitude[-1]
    gravity = MARS_G0_M_S2 * (
        MARS_RADIUS_KM / (MARS_RADIUS_KM + native_altitude[-1])
    ) ** 2
    scale_height_km = (
        KB_J_K * top_temperature_k / (mass_amu * AMU_KG * gravity) / 1000.0
    )
    result[above] = native_density[-1] * np.exp(
        -(requested_altitude[above] - native_altitude[-1]) / scale_height_km
    )
    return result


def hot_model_profile(
    native_altitude: np.ndarray,
    native_density: np.ndarray,
    requested_altitude: np.ndarray,
) -> np.ndarray:
    """Match MAMPS interpolation and return zero outside its native grid."""
    result = np.zeros(requested_altitude.shape)
    inside = (
        (requested_altitude >= native_altitude[0])
        & (requested_altitude <= native_altitude[-1])
    )
    result[inside] = np.exp(np.interp(
        requested_altitude[inside], native_altitude,
        np.log(np.maximum(native_density, 1.0e-300)),
    ))
    return result


def total_cross_section(projectile: str, target: str) -> float:
    """Reproduce MarsASPEN's internal-grid interpolation at 1000 eV."""
    table = np.loadtxt(
        CROSS_SECTIONS / f"{projectile}_{target}_cross_sections.txt",
        comments="#", skiprows=10,
    )
    columns = (2, 4, 3, 1) if projectile == "H" else (1, 2, 3, 4)
    # load_cross_sections first maps every source curve onto 768 logarithmic
    # points from 1 eV to 10 MeV. sigma_at then linearly interpolates that
    # internal grid. Repeating both stages makes this diagnostic numerically
    # consistent with Julia rather than merely close to the source table.
    internal_energy = np.logspace(0.0, 7.0, 768)
    sigma_cm2 = 0.0
    for column in columns:
        internal_sigma = np.interp(
            internal_energy, table[:, 0], table[:, column],
            left=0.0, right=0.0,
        )
        sigma_cm2 += np.interp(
            ENERGY_EV, internal_energy, internal_sigma
        )
    return sigma_cm2 * 1.0e-4


def collision_profiles(ls: int, f107: int) -> tuple[np.ndarray, np.ndarray]:
    """Return local and cumulative probabilities for H ENA and H+."""
    altitude = np.arange(80.0, 1001.0, 1.0)
    suffix = f"ls{ls:03d}_f{f107:03d}.mat"
    gitm = loadmat(ATMOSPHERE / f"gitm_{suffix}", squeeze_me=True)
    mamps = loadmat(ATMOSPHERE / f"mamps_{suffix}", squeeze_me=True)
    gitm_alt = np.asarray(gitm["alt_km"], dtype=float).squeeze()

    densities = []
    top_temperature = float(horizontal_profile(gitm, "Tn", False)[-1])
    for field, mass_amu in zip(("nCO2", "nO", "nN2"), TARGET_MASS_AMU):
        native = horizontal_profile(gitm, field, True)
        densities.append(cold_model_profile(
            gitm_alt, native, top_temperature, mass_amu, altitude
        ))

    # Transport adds MAMPS hot O to cold GITM O before evaluating n*sigma.
    hot_alt = np.asarray(mamps["alt_km"], dtype=float).squeeze()
    hot_native = horizontal_profile(mamps, "nO_hot", True)
    densities[1] += hot_model_profile(hot_alt, hot_native, altitude)

    alpha = np.zeros((2, altitude.size))
    for charge_index, (projectile, _, _) in enumerate(PROJECTILES):
        for target_index, target in enumerate(TARGETS):
            alpha[charge_index] += (
                densities[target_index]
                * total_cross_section(projectile, target)
            )

    local = 1.0 - np.exp(-alpha * 1000.0)
    # Integrate downward from 1000 km. Reversing, cumulatively summing, and
    # reversing again gives the optical depth above every altitude center.
    tau = np.cumsum((alpha[:, ::-1] * 1000.0), axis=1)[:, ::-1]
    cumulative = 1.0 - np.exp(-tau)
    return local, cumulative


def plot_grid(values_by_case: dict, cumulative: bool, output: Path) -> None:
    """Draw one 4 by 3 seasonal and solar-activity comparison."""
    altitude = np.arange(80.0, 1001.0, 1.0)
    fig, axes = plt.subplots(
        4, 3, figsize=(7.2, 8.5), sharex=True, sharey=True,
    )
    fig.subplots_adjust(
        left=0.08, right=0.99, bottom=0.06, top=0.90,
        wspace=0.11, hspace=0.18,
    )
    for row, ls in enumerate(LS_VALUES):
        for col, f107 in enumerate(F107_VALUES):
            axis = axes[row, col]
            values = values_by_case[(ls, f107)][1 if cumulative else 0]
            for index, (_, label, color) in enumerate(PROJECTILES):
                axis.plot(
                    values[index], altitude, color=color, lw=1.1, label=label
                )
            if cumulative:
                axis.set_xlim(0, 1)
            else:
                axis.set_xscale("log")
                axis.set_xlim(1.0e-10, 1)
            axis.set_ylim(80, 1000)
            axis.grid(True, which="major", color="0.90", lw=0.5)
            axis.set_title(rf"$L_s={ls}^\circ$, F10.7={f107}")
            axis.text(
                0.02, 0.97, chr(ord("a") + row * 3 + col),
                transform=axis.transAxes, ha="left", va="top",
                fontweight="bold", fontsize=8,
            )

    xlabel = (
        "Cumulative collision probability from 1000 km"
        if cumulative else "Local collision probability per 1 km"
    )
    for axis in axes[-1, :]:
        axis.set_xlabel(xlabel)
    for axis in axes[:, 0]:
        axis.set_ylabel("Altitude (km)")
    axes[0, 0].legend(ncol=2, fontsize=6, loc="lower left")
    if cumulative:
        formula = (
            r"$\tau(z)=\int_z^{1000\,\mathrm{km}}\alpha(s)\,\mathrm{d}s,"
            r"\qquad P_{\mathrm{cum}}(z)=1-\exp[-\tau(z)]$"
        )
    else:
        formula = (
            r"$\alpha(z,E)=\sum_j n_j(z)\sum_k\sigma_{j,k}(E),"
            r"\qquad P_{\mathrm{local}}(z)=1-\exp[-\alpha(z,E)\Delta s],"
            r"\quad \Delta s=1\,\mathrm{km}$"
        )
    fig.suptitle(
        rf"Fixed {ENERGY_EV:.0f} eV projectile at "
        r"lon=$0^\circ$, lat=$0^\circ$",
        fontsize=8.5,
        y=0.985,
    )
    fig.text(
        0.5, 0.955, formula,
        ha="center", va="top", fontsize=8.5,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=600, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    print(f"output={output.resolve()}")


def main() -> None:
    values = {
        (ls, f107): collision_profiles(ls, f107)
        for ls in LS_VALUES for f107 in F107_VALUES
    }
    plot_grid(
        values, False, FIGURES / "collision_probability_1000ev_local.png"
    )
    plot_grid(
        values, True, FIGURES / "collision_probability_1000ev_cumulative.png"
    )


if __name__ == "__main__":
    main()
