#!/usr/bin/env python3
'''
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
'''

""" val_report.py

"""

import os
import sys
import json
import math
import logging
import argparse
import statistics
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec



def secured_mean(a):
    try:
        result = statistics.mean(a)
    except:
        return 'NaN'

    return result

def secured_stdev(a):
    try:
        result = statistics.stdev(a)
    except:
        return 'NaN'

    return result

def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',
                        '--infile',
                        help="dam_model.json file in SI units")
    parser.add_argument('-r',
                        '--reffile',
                        help="validation_DB.json file in SI units")
    parser.add_argument('-o',
                        '--outfile',
                        help="Report.png file")
    parser.add_argument('--debug',
                        action='store_true',
                        help='Activate Debug Mode')


    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger('matplotlib').setLevel(logging.ERROR)

    args = parser.parse_args(arguments)

    logging_format = '%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s'
    if (args.debug is True):
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=logging_format)
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=logging_format)
    logging.info("Starting val_report.py")


    with open(args.infile) as m:
        model = json.load(m)

    with open(args.reffile) as r:
        ref_db = json.load(r)

    print("\ninfile =", args.infile)
    print("reffile =", args.reffile)

    damid=model["ID"]
    damname=model["Name"]
    damelev=model["Elevation"]
    Z0=model["Model"]["Z0"]
    S0=model["Model"]["S0"]
    V0=model["Model"]["V0"]
    alpha=model["Model"]["alpha"]
    beta=model["Model"]["beta"]

    print("\ndamelev =", damelev)
    print("Z0 =", Z0)
    print("S0 =", S0)
    print("V0 =", V0)
    print("alpha =", alpha)
    print("beta =", beta)

    logging.info("Model for "+damname+": S(Z) = "+format(S0, '.2F')
                 +" + "+format(alpha, '.3E')
                 +" * ( Z - "+format(Z0, '.2F')
                 +" ) ^ "+ format(beta, '.3E'))

    ref_Z0=ref_db[str(model["ID"])]["Model"]["Z0"]
    ref_S0=ref_db[str(model["ID"])]["Model"]["S0"]  # original values in m2 in DB
    ref_V0=ref_db[str(model["ID"])]["Model"]["V0"]  # original values in m3 in DB
    ref_alpha=ref_db[str(model["ID"])]["Model"]["alpha"]
    ref_beta=ref_db[str(model["ID"])]["Model"]["beta"]
    ref_Zmax=ref_db[str(model["ID"])]["Model"]["Zmax"]
    ref_Zmin=ref_db[str(model["ID"])]["Model"]["Zmin"]
    ref_Z25=ref_Zmin + 0.25 * (ref_Zmax - ref_Zmin)
    ref_Z75=ref_Zmin + 0.75 * (ref_Zmax - ref_Zmin)

    print("\nref_Z0 =", ref_Z0)
    print("ref_S0 =", ref_S0)
    print("ref_V0 =", ref_V0)
    print("ref_alpha =", ref_alpha)
    print("ref_beta =", ref_beta)

    logging.info("Reference for "+damname+": S(Z) = "+format(ref_S0, '.2F')
                 +" + "+format(ref_alpha, '.3E')
                 +" * ( Z - "+format(ref_Z0, '.2F')
                 +" ) ^ "+ format(ref_beta, '.3E'))

    print("\ndamelev  =", damelev)
    print("ref_Zmax =", ref_Zmax)

    # Figures:
    z_min = max(int(float(Z0)), int(float(ref_Z0)))
    z = range(z_min, int(float(damelev)*1.1))

    #  z_ref = range(int(float(ref_Z0)), int(float(damelev)*1.1))
    #  s_m_Zmax = S0 + alpha * math.pow((float(damelev) - Z0), beta)
    s_m_Zmax = S0 + alpha * math.pow((float(ref_Zmax) - Z0), beta)

    s_r_Zmax_m2 = ref_S0 + ref_alpha * math.pow((float(ref_Zmax) - ref_Z0), ref_beta)
    s_r_Zmin_m2 = ref_S0 + ref_alpha * math.pow((float(ref_Zmin) - ref_Z0), ref_beta)
    s_r_Z25_m2 = ref_S0 + ref_alpha * math.pow((float(ref_Z25) - ref_Z0), ref_beta)
    s_r_Z75_m2 = ref_S0 + ref_alpha * math.pow((float(ref_Z75) - ref_Z0), ref_beta)

    s = range(0, math.ceil(s_r_Zmax_m2)+10000,10000)

    Sz_model_scatter = []
    Sz_ref_scatter = []
    Szg = []
    Szh = []
    Szm = []
    Szl = []
    for h in z:
        s_m = S0 + alpha * math.pow((h - int(Z0)), beta)
        Sz_model_scatter.append(s_m)
        s_r = ref_S0 + ref_alpha * math.pow((h - int(ref_Z0)), ref_beta)
        Sz_ref_scatter.append(s_r)   #  m2

        if s_r != 0 and h >= ref_Zmin:
            Szg.append((s_r-s_m)/(s_r))
            if h >= ref_Z75:
                Szh.append((s_r-s_m)/(s_r))
            elif h >= ref_Z25:
                Szm.append((s_r-s_m)/(s_r))
            else:
                Szl.append((s_r-s_m)/(s_r))


    print("\n== s_m =",s_m)
    print("== s_r =",s_r)

    print("\n== Zmax model =",s_m_Zmax)
    print("== Zmax ref   =", s_r_Zmax_m2)

    Vs_model_scatter = []
    Vs_ref_scatter = []

    v_m_Zmax = V0 + math.pow(((s_r_Zmax_m2-S0)/alpha), 1/beta) * (S0 + ((s_r_Zmax_m2-S0) / (beta + 1.)))
    v_r_Zmax = ref_V0 + math.pow(((s_r_Zmax_m2-ref_S0)/ref_alpha), 1/ref_beta) * (ref_S0 + ((s_r_Zmax_m2-ref_S0) / (ref_beta + 1.)))

    Tx_model_scatter = []
    Tx_ref_scatter = []
    Vsg = []
    Vsh = []
    Vsm = []
    Vsl = []
    Tsg = []
    Tsh = []
    Tsm = []
    Tsl = []
    #  s_m_norm = []
    #  s_r_norm = []
    for a in s:
        v_m = V0 + math.pow(((a-S0)/alpha), 1/beta) * (S0 + ((a-S0) / (beta + 1.)))
        Vs_model_scatter.append(v_m)
        Tx_model_scatter.append(v_m/v_m_Zmax)
        #  s_m_norm.append(a/s_m_max)

        #  v_r = ref_V0 + math.pow(((a*10000-ref_S0*10000)/ref_alpha), 1/ref_beta) * (ref_S0*10000 + ((a*10000-ref_S0*10000) / (ref_beta + 1.)))
        v_r = ref_V0 + math.pow(((a-ref_S0)/ref_alpha), 1/ref_beta) * (ref_S0 + ((a-ref_S0) / (ref_beta + 1.)))
        Vs_ref_scatter.append(v_r)
        Tx_ref_scatter.append(v_r/(v_r_Zmax))

        #  s_r_norm.append(a/(s_r_max*10000))

        if v_r != 0 and a >= s_r_Zmin_m2:
            Vsg.append((v_r-v_m)/(v_r))
            Tsg.append(((v_r/v_r_Zmax)-(v_m/v_m_Zmax))/(v_r/v_r_Zmax))
            if a >= s_r_Z75_m2:
                Vsh.append((v_r-v_m)/(v_r))
                Tsh.append(((v_r/v_r_Zmax)-(v_m/v_m_Zmax))/(v_r/v_r_Zmax))
            if a >= s_r_Z25_m2:
                Vsm.append((v_r-v_m)/(v_r))
                Tsm.append(((v_r/v_r_Zmax)-(v_m/v_m_Zmax))/(v_r/v_r_Zmax))
            else:
                Vsl.append((v_r-v_m)/(v_r))
                Tsl.append(((v_r/v_r_Zmax)-(v_m/v_m_Zmax))/(v_r/v_r_Zmax))

    print("\n== v_m_Zmax =",v_m_Zmax)
    print("== v_r_Zmax =",v_r_Zmax)

    ticks_m2 = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/10000.))
    ticks_m3 = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000000.))

    fig = plt.figure(dpi=300)
    #  fig.subplots_adjust(hspace=0)
    gs = fig.add_gridspec( 2, 1, height_ratios=[1, 1])

    Sz = plt.subplot(gs[0])
    Sz.plot(Sz_ref_scatter, Sz_ref_scatter,
            'r-')
    Sz.scatter(Sz_ref_scatter, Sz_model_scatter,
               label='S(z)_model = f(S(z)_reference)')
    Sz.grid(b=True, which='major', linestyle='-')
    Sz.grid(b=False, which='minor', linestyle='--')
    Sz.set_xlabel('S(z)_reference (ha)', fontsize=6)
    Sz.set_ylabel('S(z)_model (ha)', fontsize=6)
    Sz.tick_params(axis='both', which='major', labelsize=6)
    # Trick to display in Ha
    Sz.xaxis.set_major_formatter(ticks_m2)
    Sz.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    Sz_comp = plt.subplot(gs[1])
    Sz_comp.plot(z, Sz_ref_scatter,
            'r-',
            label='S(z)_reference')
    Sz_comp.plot(z, Sz_model_scatter,
            'g-',
            label='S(z)_model')
    Sz_comp.axvline(x=float(ref_Zmax),
                    ls=':', lw=2,
                    color='red',
                    label='Reference Max Dam Z')
    Sz_comp.axvline(x=float(damelev),
                    ls=':', lw=2,
                    color='green',
                    label='Model Max Dam Z')
    Sz_comp.grid(b=True, which='major', linestyle='-')
    Sz_comp.grid(b=False, which='minor', linestyle='--')
    Sz_comp.set_xlabel('z (m)', fontsize=6)
    Sz_comp.set_ylabel('S(z) (ha)', fontsize=6)
    Sz_comp.tick_params(axis='both', which='major', labelsize=6)
    # Trick to display in Ha
    Sz_comp.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    plt.suptitle(damname + ": S(z)",
                 fontsize=10)
    plt.savefig(os.path.splitext(args.outfile)[0]+"_Sz.png", dpi=300)


    fig = plt.figure(dpi=300)
    #  fig.subplots_adjust(hspace=0)
    gv = fig.add_gridspec( 2, 1, height_ratios=[1, 1])

    Vs = plt.subplot(gv[0])
    Vs.plot(Vs_ref_scatter, Vs_ref_scatter,
            'r-')
    Vs.scatter(Vs_ref_scatter, Vs_model_scatter,
               label='V(S)_model = f(V(S)_reference)')
    Vs.grid(b=True, which='major', linestyle='-')
    Vs.grid(b=False, which='minor', linestyle='--')
    Vs.set_xlabel('V(S)_reference (hm3)', fontsize=6)
    Vs.set_ylabel('V(S)_model (hm3)', fontsize=6)
    Vs.tick_params(axis='both', which='major', labelsize=6)
    # Trick to display in hm3
    Vs.xaxis.set_major_formatter(ticks_m3)
    Vs.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    Vs_comp = plt.subplot(gs[1])
    Vs_comp.plot(s, Vs_ref_scatter,
            'r-',
            label='V(S)_reference')
    Vs_comp.plot(s, Vs_model_scatter,
            'g-',
            label='V(S)_model')
    Vs_comp.axvline(x=s_r_Zmax_m2,
                    ls=':', lw=2,
                    color='red',
                    label='S(Zmax_ref)')
    Vs_comp.grid(b=True, which='major', linestyle='-')
    Vs_comp.grid(b=False, which='minor', linestyle='--')
    Vs_comp.set_xlabel('S (ha)', fontsize=6)
    Vs_comp.set_ylabel('V(S) (hm3)', fontsize=6)
    Vs_comp.tick_params(axis='both', which='major', labelsize=6)
    # Trick to display in Ha
    Vs_comp.xaxis.set_major_formatter(ticks_m2)
    Vs_comp.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    plt.suptitle(damname + ": V(S)",
                 fontsize=10)
    plt.savefig(os.path.splitext(args.outfile)[0]+"_Vs.png", dpi=300)

    fig = plt.figure(dpi=300)
    #  fig.subplots_adjust(hspace=0)
    gt = fig.add_gridspec( 2, 1, height_ratios=[1, 1])

    Tx = plt.subplot(gt[0])
    Tx.plot(Tx_ref_scatter, Tx_ref_scatter,
            'r-')
    Tx.scatter(Tx_ref_scatter, Tx_model_scatter,
               label='Vr(S)_model = f(Vr(S)_reference)')
    Tx.grid(b=True, which='major', linestyle='-')
    Tx.grid(b=False, which='minor', linestyle='--')
    Tx.set_xlabel('Vr(S)_reference (%)', fontsize=6)
    Tx.set_ylabel('Vr(S)_model (%)', fontsize=6)
    Tx.tick_params(axis='both', which='major', labelsize=6)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    Tx_comp = plt.subplot(gt[1])
    Tx_comp.plot(s, Tx_ref_scatter,
            'r-',
            label='Vr(S)_reference')
    Tx_comp.plot(s, Tx_model_scatter,
            'g-',
            label='Vr(S)_model')
    Tx_comp.axvline(x=s_r_Zmax_m2,
                    ls=':', lw=2,
                    color='red',
                    label='S(Zmax_ref)')
    Tx_comp.grid(b=True, which='major', linestyle='-')
    Tx_comp.grid(b=False, which='minor', linestyle='--')
    Tx_comp.set_xlabel('S (ha)', fontsize=6)
    Tx_comp.set_ylabel('Vr(S) (%)', fontsize=6)
    Tx_comp.tick_params(axis='both', which='major', labelsize=6)
    # Trick to display in Ha
    Tx_comp.xaxis.set_major_formatter(ticks_m2)
    Tx_comp.set_ylim(-0.1, 1.1)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    plt.suptitle(damname + ": Volume Rate",
                 fontsize=10)
    plt.savefig(os.path.splitext(args.outfile)[0]+"_VolumeRate.png", dpi=300)


    results_json = {
        'ID': damid,
        'Name': damname,
        'Zmin': ref_Zmin,
        'Z_25': ref_Z25,
        'Z_75': ref_Z75,
        'Zmax': ref_Zmax,
        'Smin': s_r_Zmin_m2,
        'S_25': s_r_Z25_m2,
        'S_75': s_r_Z75_m2,
        'Smax': s_r_Zmax_m2,
        'S(z)_quality': {
            'glob' : {
                'mean' : secured_mean(Szg),
                'stdev' : secured_stdev(Szg)
            },
            'high' : {
                'mean' : secured_mean(Szh),
                'stdev' : secured_stdev(Szh)
            },
            'mid' : {
                'mean' : secured_mean(Szm),
                'stdev' : secured_stdev(Szm)
            },
            'low' : {
                'mean' : secured_mean(Szl),
                'stdev' : secured_stdev(Szl)
            }
        },
        'V(S)_quality': {
            'glob' : {
                'mean' : secured_mean(Vsg),
                'stdev' : secured_stdev(Vsg)
            },
            'high' : {
                'mean' : secured_mean(Vsh),
                'stdev' : secured_stdev(Vsh)
            },
            'mid' : {
                'mean' : secured_mean(Vsm),
                'stdev' : secured_stdev(Vsm)
            },
            'low' : {
                'mean' : secured_mean(Vsl),
                'stdev' : secured_stdev(Vsl)
            }
        },
        'Vr(S)_quality': {
            'glob' : {
                'mean' : secured_mean(Tsg),
                'stdev' : secured_stdev(Tsg)
            },
            'high' : {
                'mean' : secured_mean(Tsh),
                'stdev' : secured_stdev(Tsh)
            },
            'mid' : {
                'mean' : secured_mean(Tsm),
                'stdev' : secured_stdev(Tsm)
            },
            'low' : {
                'mean' : secured_mean(Tsl),
                'stdev' : secured_stdev(Tsl)
            }
        },
        'Dam_bottom_estimation' : {
            'Z0_mod' : Z0,
            'Z0_ref' : ref_Z0,
            'Dam_bottom_error' : abs(Z0 - ref_Z0)
        }
    }

    with open(os.path.splitext(args.outfile)[0]+".json", "w") as write_file:
        json.dump(results_json, write_file, indent=4)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
