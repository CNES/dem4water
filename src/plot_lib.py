#!/usr/bin/env python3
"""
:author: Benjamin Tardy <benjamin.tardy@csgroup.eu>
:organization: CS Group
:copyright: 2022 CNES. All rights reserved.
:license: see LICENSE file
:created: 2022
"""
import matplotlib.pyplot as plt


def plot_slope(all_zi, all_szi, model_zi, model_szi, model_mae, damelev, abs_zi, abs_szi, abs_mae):
    """"""
    alt = range(int(model_zi[0])+1, int(model_zi[-1]))
    ms_fig = plt.figure(dpi=300)
    ms_fig.subplots_adjust(hspace=0)
    ms_gs = ms_fig.add_gridspec(3, 1, height_ratios=[1, 1, 1])

    maeax = plt.subplot(ms_gs[0])
    maeax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    maeax.plot(model_zi, model_mae, marker=".", color="purple", linestyle="dashed", label="MAE")
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
        median(model_Zi),
        best_P[0],
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
    ds = np.diff(l_slope) / np.diff(l_z)
    dz = (np.array(l_z)[:-1] + np.array(l_z)[1:]) / 2
    dslpax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    dslpax.plot(dz, ds, marker=".", color="blue", label="dSlope")
    dslpax.grid(b=True, which="major", linestyle="-")
    dslpax.grid(b=False, which="minor", linestyle="--")
    dslpax.set(
        xlabel="Virtual Water Surface Elevation (m)", ylabel="Local Slope Derivative"
    )
    dslpax.set_xlim(alt[0], alt[-1] + 10)
    # Trick to display in Ha
    plt.minorticks_on()

    # if args.debug is True:
    plt.savefig(os.path.splitext(args.outfile)[0] + "_slope.png", dpi=300)

def plot_model_combo():
    """"""
    fig = plt.figure(dpi=300)
    fig.subplots_adjust(hspace=0)
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1])

    maeax0 = plt.subplot(gs[0])
    maeax0.axvline(x=float(damelev), ls=":", lw=2, color="teal", label="Dam Elevation")
    maeax0.plot(z, abs_sz, "g--", label="S(z) abs_model")
    maeax0.plot(z, sz, "r--", label="S(z)")
    maeax0.plot(Zi[0], S_Zi[0], "ro")
    maeax0.scatter(Zi, S_Zi, marker=".", label="S(Z_i)")
    maeax0.plot(
        Zi[best_i : best_i + args.winsize],
        S_Zi[best_i : best_i + args.winsize],
        marker="o",
        ls=":",
        lw=1,
        color="yellow",
        markerfacecolor="blue",
        label="Selected S(Z_i)",
    )
    maeax0.plot(
        median(Zi[best_i : best_i + args.winsize]),
        median(S_Zi[best_i : best_i + args.winsize]),
        color="#ff7f0e",
        marker="p",
        markerfacecolor="#ff7f0e",
    )
    maeax0.grid(b=True, which="major", linestyle="-")
    maeax0.grid(b=False, which="minor", linestyle="--")
    maeax0.set(
        xlabel="Virtual Water Surface Elevation (m)",
        ylabel="Virtual Water Surface (ha)",
    )
    maeax0.label_outer()
    #  maeax0.set_xlim(z[0]-10, z[-1]+10)
    if data_shortage is False:
        maeax0.set_xlim(z[0] - 5, median(Zi[best_i : best_i + args.winsize]) + 30)
        maeax0.set_ylim(-10, S_Zi[best_i + args.winsize] * 1.1)
    # Trick to display in Ha
    maeax0.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")

    maeax1 = plt.subplot(gs[1])
    maeax1.axvline(x=float(damelev), ls=":", lw=2, color="teal", label="Dam Elevation")
    maeax1.plot(l_z, l_mae, color="purple", marker=".", linestyle="dashed", label="MAE")
    maeax1.plot(
        median(Zi[best_i : best_i + args.winsize]),
        best,
        "p",
        color="#ff7f0e",
        label="Selected MAE",
    )
    maeax1.plot(
        median(Zi[abs_i : abs_i + args.winsize]),
        abs_mae,
        "+",
        color="black",
        label="Absolute min(MAE)",
    )
    maeax1.grid(b=True, which="major", linestyle="-")
    maeax1.grid(b=True, which="minor", linestyle="--")
    maeax1.set(xlabel="Virtual Water Surface Elevation (m)", ylabel="Local MAE (ha)")
    #  maeax1.set_xlim(z[0]-10, z[-1]+10)
    if data_shortage is False:
        maeax1.set_xlim(z[0] - 5, median(Zi[best_i : best_i + args.winsize]) + 30)
    maeax1.set_yscale("log")
    # Trick to display in Ha
    maeax1.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="lower left")

    plt.suptitle(
        damname
        + ": S(Z) = "
        + format(S_Zi[0], ".2F")
        + " + "
        + format(best_alpha, ".3E")
        + " * ( Z - "
        + format(Zi[0], ".2F")
        + " ) ^ "
        + format(best_beta, ".3E"),
        fontsize=10,
    )
    plt.savefig(os.path.splitext(args.outfile)[0] + "_combo.png", dpi=300)

def plot_model(all_szi, all_zi, model_szi,model_zi, alpha, beta, outfile):
    """"""
    alt = range(int(model_zi[0])+1, int(model_zi[-1]))
    mod_sz = []
    for curr_alt in alt:
        mod_sz.append(model_szi[0] + alpha * math.pow((curr_alt - model_zi[0]), beta))
    
    fig, ax = plt.subplots()
    ax.plot(alt, mod_sz, "r--", label="S(z)")
    ax.plot(model_zi[0], model_szi[0], "ro")
    ax.scatter(model_zi[1:], model_szi[1:], marker=".", label="S(Zi)")
    ax.scatter(
        median(model_zi),
        median(model_szi),
        color="#ff7f0e",
        marker="p",
        label="Selected S(Zi)",
    )
    #  ax.plot(z, reg_abs(z), 'b-')
    #  ax.plot(z, reg_best(z), 'g-')
    ax.grid(visible=True, which="major", linestyle="-")
    ax.grid(visible=False, which="minor", linestyle="--")
    ax.set(
        xlabel="Virtual Water Surface Elevation (m)",
        ylabel="Virtual Water Surface (ha)",
    )
    plt.title(
        damname
        + ": S(Z) = "
        + format(model_szi[0], ".2F")
        + " + "
        + format(alpha, ".3E")
        + " * ( Z - "
        + format(model_zi[0], ".2F")
        + " ) ^ "
        + format(beta, ".3E"),
        fontsize=10,
    )
    ax.set_xlim(alt[0] - 10, alt[-1] + 10)
    # Trick to display in Ha
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 10000.0))
    ax.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")
    plt.savefig(outfile, dpi=300)

def plot_vs():
    S_dam = S_Zi[0] + best_alpha * math.pow((float(damelev) - Zi[0]), best_beta)
    s = range(0, math.ceil(S_dam) + 10000)
    v0 = 0.0
    vs = []
    for a in s:
        v = v0 + math.pow(((a - S_Zi[0]) / best_alpha), 1 / best_beta) * (
            S_Zi[0] + ((a - S_Zi[0]) / (best_beta + 1.0))
        )
        vs.append(v)

    # V(S) Moldel Plot
    figv, axv = plt.subplots()
    axv.axvline(x=S_dam, ls=":", lw=2, color="teal", label="Max Dam Surface")
    axv.plot(s, vs, "r-", label="V(S)")
    axv.grid(b=True, which="major", linestyle="-")
    axv.grid(b=False, which="minor", linestyle="--")
    axv.set(
        xlabel="Estimated Water Surface (ha)", ylabel="Modelized Water Volume (hm3)"
    )
    plt.title(
        damname + " : V=f(S) ",
        #  +": S(Z) = "+format(S_Zi[0], '.2F')+" + "+format(best_alpha, '.3E')
        #  +" * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'),
        fontsize=10,
    )
    # Trick to display in Ha
    axv.xaxis.set_major_formatter(ticks_m2)
    ticks_m3 = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000.0))
    axv.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={"size": 6}, loc="upper left")
    plt.savefig(os.path.splitext(args.outfile)[0] + "_VS.png", dpi=300)
