#!/usr/bin/env python3
"""Graph module."""
import math
from statistics import median

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker


def plot_szi_points(r_elev, r_area, pdb_elev, damname, outfile):
    """Plot the szi using area and altitudes."""
    fig, ax1 = plt.subplots()
    # Trick to display in Ha
    ticks_y = ticker.FuncFormatter(lambda x, pos: f"{x/10000.0:g}")
    ax1.yaxis.set_major_formatter(ticks_y)
    ax1.plot(
        r_elev[:-1],
        r_area[:-1],
        color="#ff7f0e",
        marker="o",
        linestyle="dashed",
        markerfacecolor="blue",
    )
    ax1.set_ylim(bottom=0.0)
    ax1.set(
        xlabel="Virtual Water Surface Elevation (m)",
        ylabel="Virtual Water Surface (ha)",
    )
    ax1.plot(pdb_elev, 0.0, "ro")
    ax1.grid(visible=True, which="major", linestyle="-")
    ax1.grid(visible=False, which="minor", linestyle="--")
    plt.minorticks_on()
    plt.title(damname + ": S(Z_i)")
    fig.savefig(outfile)


def plot_slope(
    model_zi,
    model_mae,
    damelev,
    abs_zi,
    abs_mae,
    best,
    damname,
    l_z,
    l_slope,
    best_p,
    outfile,
):
    """Plot slope graph."""
    alt = range(int(model_zi[0]) + 1, int(model_zi[-1]))
    ms_fig = plt.figure(dpi=300)
    ms_fig.subplots_adjust(hspace=0)
    ms_gs = ms_fig.add_gridspec(3, 1, height_ratios=[1, 1, 1])

    maeax = plt.subplot(ms_gs[0])
    maeax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    maeax.plot(
        l_z, model_mae, marker=".", color="purple", linestyle="dashed", label="MAE"
    )
    maeax.plot(
        median(model_zi),
        best,
        "p",
        color="#ff7f0e",
        label="Selected min(MAE)",
    )
    maeax.plot(
        median(abs_zi),
        abs_mae,
        "+",
        color="black",
        label="Absolute min(MAE)",
    )
    maeax.grid(visible=True, which="major", linestyle="-")
    maeax.grid(visible=False, which="minor", linestyle="--")
    maeax.set(xlabel="Virtual Water Surface Elevation (m)", ylabel="Local MAE (ha)")
    maeax.set_xlim(alt[0], alt[-1] + 10)
    maeax.set_yscale("log")
    maeax.label_outer()
    # Trick to display in Ha
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: f"{x / 10000.0:g}")
    maeax.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.title(damname + ": Local Maximum Absolute Error")
    plt.legend(prop={"size": 6})

    # Plot slope
    #  slpfig, slpax = plt.subplots()
    slpax = plt.subplot(ms_gs[1])
    slpax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    slpax.plot(
        l_z, l_slope, marker=".", color="purple", linestyle="dashed", label="Slope"
    )
    slpax.plot(
        median(model_zi),
        best_p[0],
        "p",
        color="#ff7f0e",
        label="Selected min(MAE)",
    )
    slpax.grid(b=True, which="major", linestyle="-")
    slpax.grid(b=False, which="minor", linestyle="--")
    slpax.set(xlabel="Virtual Water Surface Elevation (m)", ylabel="Local Slope")
    slpax.set_xlim(alt[0], alt[-1] + 10)
    # Trick to display in Ha
    plt.minorticks_on()
    dslpax = plt.subplot(ms_gs[2])
    der_s = np.diff(l_slope) / np.diff(l_z)
    der_z = (np.array(l_z)[:-1] + np.array(l_z)[1:]) / 2
    dslpax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    dslpax.plot(der_z, der_s, marker=".", color="blue", label="dSlope")
    dslpax.grid(visible=True, which="major", linestyle="-")
    dslpax.grid(visible=False, which="minor", linestyle="--")
    dslpax.set(
        xlabel="Virtual Water Surface Elevation (m)", ylabel="Local Slope Derivative"
    )
    dslpax.set_xlim(alt[0], alt[-1] + 10)
    # Trick to display in Ha
    plt.minorticks_on()

    # if args.debug is True:
    plt.savefig(outfile, dpi=300)


def plot_model_combo(
    all_zi,
    all_szi,
    model_zi,
    model_szi,
    damelev,
    data_shortage,
    alpha,
    beta,
    damname,
    alt_sz,
    abs_sz,
    mod_sz,
    l_z,
    l_mae,
    abs_mae,
    best,
    outfile,
):
    """Plot combo graph."""
    fig = plt.figure(dpi=300)
    fig.subplots_adjust(hspace=0)
    g_s = fig.add_gridspec(2, 1, height_ratios=[4, 1])

    alt = range(int(all_zi[0]) + 1, int(all_zi[-1]))
    print(alt)
    maeax0 = plt.subplot(g_s[0])
    maeax0.axvline(x=float(damelev), ls=":", lw=2, color="teal", label="Dam Elevation")
    maeax0.plot(alt_sz, abs_sz, "g--", label="S(z) abs_model")
    maeax0.plot(alt_sz, mod_sz, "r--", label="S(z)")
    maeax0.plot(all_zi[0], all_szi[0], "ro")
    maeax0.scatter(all_zi, all_szi, marker=".", label="S(Z_i)")
    maeax0.plot(
        model_zi,
        model_szi,
        marker="o",
        ls=":",
        lw=1,
        color="yellow",
        markerfacecolor="blue",
        label="Selected S(Z_i)",
    )
    maeax0.plot(
        median(model_zi),
        median(model_szi),
        color="#ff7f0e",
        marker="p",
        markerfacecolor="#ff7f0e",
    )
    maeax0.grid(visible=True, which="major", linestyle="-")
    maeax0.grid(visible=False, which="minor", linestyle="--")
    maeax0.set(
        xlabel="Virtual Water Surface Elevation (m)",
        ylabel="Virtual Water Surface (ha)",
    )
    maeax0.label_outer()
    #  maeax0.set_xlim(z[0]-10, z[-1]+10)
    maeax0.set_xlim(alt[0] - 10, alt[-1] + 10)
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: f"{x/10000.0:g}")
    if data_shortage is False:
        maeax0.set_xlim(alt[0] - 5, median(model_zi) + 30)
        maeax0.set_ylim(-10, model_szi[-1] * 1.1)
    # Trick to display in Ha
    maeax0.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    maeax1 = plt.subplot(g_s[1])
    maeax1.axvline(x=float(damelev), ls=":", lw=2, color="teal", label="Dam Elevation")
    maeax1.plot(l_z, l_mae, color="purple", marker=".", linestyle="dashed", label="MAE")
    maeax1.plot(
        median(model_zi),
        best,
        "p",
        color="#ff7f0e",
        label="Selected MAE",
    )
    maeax1.plot(
        median(model_zi),
        abs_mae,
        "+",
        color="black",
        label="Absolute min(MAE)",
    )
    maeax1.grid(visible=True, which="major", linestyle="-")
    maeax1.grid(visible=True, which="minor", linestyle="--")
    maeax1.set(xlabel="Virtual Water Surface Elevation (m)", ylabel="Local MAE (ha)")
    maeax1.set_xlim(alt[0] - 10, alt[-1] + 10)
    if data_shortage is False:
        maeax1.set_xlim(alt[0] - 5, median(model_zi) + 30)
    maeax1.set_yscale("log")
    # Trick to display in Ha
    maeax1.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="lower left")

    plt.suptitle(
        f"{damname}: S(Z) = {all_szi[0]:.2F} + {alpha:.3E}"
        f" * ( Z - {all_zi[0]:.2F}  ) ^ {beta:.3E}",
        fontsize=10,
    )
    plt.savefig(outfile, dpi=300)


def plot_model(model_szi, model_zi, z_0, sz_0, alpha, beta, damname, outfile):
    """Plot model."""
    alt = range(int(model_zi[0]) + 1, int(model_zi[-1]))
    mod_sz = []
    for curr_alt in alt:
        mod_sz.append(sz_0 + alpha * math.pow((curr_alt - z_0), beta))

    _, axe = plt.subplots()
    axe.plot(alt, mod_sz, "r--", label="S(z)")
    axe.plot(z_0, sz_0, "ro")
    axe.scatter(model_zi[1:], model_szi[1:], marker=".", label="S(Zi)")
    axe.scatter(
        median(model_zi),
        median(model_szi),
        color="#ff7f0e",
        marker="p",
        label="Selected S(Zi)",
    )
    #  ax.plot(z, reg_abs(z), 'b-')
    #  ax.plot(z, reg_best(z), 'g-')
    axe.grid(visible=True, which="major", linestyle="-")
    axe.grid(visible=False, which="minor", linestyle="--")
    axe.set(
        xlabel="Virtual Water Surface Elevation (m)",
        ylabel="Virtual Water Surface (ha)",
    )
    plt.title(
        damname
        + ": S(Z) = "
        + format(sz_0, ".2F")
        + " + "
        + format(alpha, ".3E")
        + " * ( Z - "
        + format(z_0, ".2F")
        + " ) ^ "
        + format(beta, ".3E"),
        fontsize=10,
    )
    axe.set_xlim(z_0 - 10, alt[-1] + 10)
    # Trick to display in Ha
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: f"{x / 10000.0:g}")
    axe.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")
    plt.savefig(outfile, dpi=300)


def plot_vs(all_zi, all_szi, damelev, alpha, beta, damname, outfile):
    """Plot V(S)."""
    s_dam = all_szi[0] + alpha * math.pow((float(damelev) - all_zi[0]), beta)
    surf = range(0, math.ceil(s_dam) + 10000)
    vol0 = 0.0
    vol_s = []
    for area in surf:
        volume = vol0 + math.pow(((area - all_szi[0]) / alpha), 1 / beta) * (
            all_szi[0] + ((area - all_szi[0]) / (beta + 1.0))
        )
        vol_s.append(volume)

    # V(S) Moldel Plot
    _, axv = plt.subplots()
    axv.axvline(x=s_dam, ls=":", lw=2, color="teal", label="Max Dam Surface")
    axv.plot(surf, vol_s, "r-", label="V(S)")
    axv.grid(visible=True, which="major", linestyle="-")
    axv.grid(visible=False, which="minor", linestyle="--")
    axv.set(
        xlabel="Estimated Water Surface (ha)", ylabel="Modelized Water Volume (hm3)"
    )
    plt.title(
        f"{damname} : V=f(S) ",
        fontsize=10,
    )
    #  +": S(Z) = "+format(S_Zi[0], '.2F')+" + "+format(best_alpha, '.3E')
    #  +" * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'),
    # Trick to display in Ha
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: f"{x/10000.0:g}")
    axv.xaxis.set_major_formatter(ticks_m2)
    ticks_m3 = ticker.FuncFormatter(lambda x, pos: f"{x / 1000000.0:g}")
    axv.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")
    plt.savefig(outfile, dpi=300)


def plot_report_sz(
    alt, sz_ref_scatter, sz_model_scatter, ref_zmax, damelev, damname, outfile
):
    """Plot Sz ref and model."""
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: f"{x/10000.0:g}")
    fig = plt.figure(dpi=300)
    #  fig.subplots_adjust(hspace=0)
    g_s = fig.add_gridspec(2, 1, height_ratios=[1, 1])

    sz1 = plt.subplot(g_s[0])
    sz1.plot(sz_ref_scatter, sz_ref_scatter, "r-")
    sz1.scatter(
        sz_ref_scatter, sz_model_scatter, label="S(z)_model = f(S(z)_reference)"
    )
    sz1.grid(b=True, which="major", linestyle="-")
    sz1.grid(b=False, which="minor", linestyle="--")
    sz1.set_xlabel("S(z)_reference (ha)", fontsize=6)
    sz1.set_ylabel("S(z)_model (ha)", fontsize=6)
    sz1.tick_params(axis="both", which="major", labelsize=6)
    # Trick to display in Ha
    sz1.xaxis.set_major_formatter(ticks_m2)
    sz1.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    sz_comp = plt.subplot(g_s[1])
    sz_comp.plot(alt, sz_ref_scatter, "r-", label="S(z)_reference")
    sz_comp.plot(alt, sz_model_scatter, "g-", label="S(z)_model")
    sz_comp.axvline(
        x=float(ref_zmax), ls=":", lw=2, color="red", label="Reference Max Dam Z"
    )
    sz_comp.axvline(
        x=float(damelev), ls=":", lw=2, color="green", label="Model Max Dam Z"
    )
    sz_comp.grid(b=True, which="major", linestyle="-")
    sz_comp.grid(b=False, which="minor", linestyle="--")
    sz_comp.set_xlabel("z (m)", fontsize=6)
    sz_comp.set_ylabel("S(z) (ha)", fontsize=6)
    sz_comp.tick_params(axis="both", which="major", labelsize=6)
    # Trick to display in Ha
    sz_comp.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    plt.suptitle(damname + ": S(z)", fontsize=10)
    plt.savefig(outfile, dpi=300)


def plot_report_vs(
    surf, vs_ref_scatter, vs_model_scatter, s_r_zmax_m2, damname, outfile
):
    """Plot VS report."""
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: f"{x/10000.0:g}")
    ticks_m3 = ticker.FuncFormatter(lambda x, pos: f"{x / 1000000.0:g}")
    fig = plt.figure(dpi=300)
    #  fig.subplots_adjust(hspace=0)
    g_v = fig.add_gridspec(2, 1, height_ratios=[1, 1])

    vs1 = plt.subplot(g_v[0])
    vs1.plot(vs_ref_scatter, vs_ref_scatter, "r-")
    vs1.scatter(
        vs_ref_scatter, vs_model_scatter, label="V(S)_model = f(V(S)_reference)"
    )
    vs1.grid(b=True, which="major", linestyle="-")
    vs1.grid(b=False, which="minor", linestyle="--")
    vs1.set_xlabel("V(S)_reference (hm3)", fontsize=6)
    vs1.set_ylabel("V(S)_model (hm3)", fontsize=6)
    vs1.tick_params(axis="both", which="major", labelsize=6)
    # Trick to display in hm3
    vs1.xaxis.set_major_formatter(ticks_m3)
    vs1.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    vs_comp = plt.subplot(g_v[1])
    vs_comp.plot(surf, vs_ref_scatter, "r-", label="V(S)_reference")
    vs_comp.plot(surf, vs_model_scatter, "g-", label="V(S)_model")
    vs_comp.axvline(x=s_r_zmax_m2, ls=":", lw=2, color="red", label="S(Zmax_ref)")
    vs_comp.grid(visible=True, which="major", linestyle="-")
    vs_comp.grid(visible=False, which="minor", linestyle="--")
    vs_comp.set_xlabel("S (ha)", fontsize=6)
    vs_comp.set_ylabel("V(S) (hm3)", fontsize=6)
    vs_comp.tick_params(axis="both", which="major", labelsize=6)
    # Trick to display in Ha
    vs_comp.xaxis.set_major_formatter(ticks_m2)
    vs_comp.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    plt.suptitle(damname + ": V(S)", fontsize=10)
    plt.savefig(outfile, dpi=300)


def plot_report_volume_rate(
    surf, tx_ref_scatter, tx_model_scatter, s_r_zmax_m2, damname, outfile
):
    """Plot Volume Rate report."""
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: f"{x/10000.0:g}")
    fig = plt.figure(dpi=300)
    #  fig.subplots_adjust(hspace=0)
    gt = fig.add_gridspec(2, 1, height_ratios=[1, 1])

    tx1 = plt.subplot(gt[0])
    tx1.plot(tx_ref_scatter, tx_ref_scatter, "r-")
    tx1.scatter(
        tx_ref_scatter, tx_model_scatter, label="Vr(S)_model = f(Vr(S)_reference)"
    )
    tx1.grid(b=True, which="major", linestyle="-")
    tx1.grid(b=False, which="minor", linestyle="--")
    tx1.set_xlabel("Vr(S)_reference (%)", fontsize=6)
    tx1.set_ylabel("Vr(S)_model (%)", fontsize=6)
    tx1.tick_params(axis="both", which="major", labelsize=6)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    tx_comp = plt.subplot(gt[1])
    tx_comp.plot(surf, tx_ref_scatter, "r-", label="Vr(S)_reference")
    tx_comp.plot(surf, tx_model_scatter, "g-", label="Vr(S)_model")
    tx_comp.axvline(x=s_r_zmax_m2, ls=":", lw=2, color="red", label="S(Zmax_ref)")
    tx_comp.grid(b=True, which="major", linestyle="-")
    tx_comp.grid(b=False, which="minor", linestyle="--")
    tx_comp.set_xlabel("S (ha)", fontsize=6)
    tx_comp.set_ylabel("Vr(S) (%)", fontsize=6)
    tx_comp.tick_params(axis="both", which="major", labelsize=6)
    # Trick to display in Ha
    tx_comp.xaxis.set_major_formatter(ticks_m2)
    tx_comp.set_ylim(-0.1, 1.1)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    plt.suptitle(damname + ": Volume Rate", fontsize=10)
    plt.savefig(outfile, dpi=300)
