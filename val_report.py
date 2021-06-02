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
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',
                        '--infile',
                        help="dam_model.json file")
    parser.add_argument('-r',
                        '--reffile',
                        help="validation_DB.json file")
    parser.add_argument('-o',
                        '--outfile',
                        help="Repport.png file")
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

    damname=model["Name"]
    damelev=model["Elevation"]
    Z0=model["Model"]["Z0"]
    S0=model["Model"]["S0"]
    V0=model["Model"]["V0"]
    alpha=model["Model"]["alpha"]
    beta=model["Model"]["beta"]

    logging.info("Model for "+damname+": S(Z) = "+format(S0, '.2F')
                 +" + "+format(alpha, '.3E')
                 +" * ( Z - "+format(Z0, '.2F')
                 +" ) ^ "+ format(beta, '.3E'))

    ref_Z0=ref_db[str(model["ID"])]["Model"]["Z0"]
    ref_S0=ref_db[str(model["ID"])]["Model"]["S0"]
    ref_V0=ref_db[str(model["ID"])]["Model"]["V0"]
    ref_alpha=ref_db[str(model["ID"])]["Model"]["alpha"]
    ref_beta=ref_db[str(model["ID"])]["Model"]["beta"]
    ref_Zmax=ref_db[str(model["ID"])]["Model"]["Zmax"]

    logging.info("Reference for "+damname+": S(Z) = "+format(ref_S0, '.2F')
                 +" + "+format(ref_alpha, '.3E')
                 +" * ( Z - "+format(ref_Z0, '.2F')
                 +" ) ^ "+ format(ref_beta, '.3E'))

    # Figures:
    z_min = max(int(float(Z0)), int(float(ref_Z0)))
    z = range(z_min, int(float(damelev)*1.1))
    z_mod = range(int(float(Z0)), int(float(damelev)*1.1))
    z_ref = range(int(float(ref_Z0)), int(float(damelev)*1.1))
    S_dam = S0 + alpha * math.pow((float(damelev) - Z0), beta)
    s = range(0, math.ceil(S_dam)+10000,10000)

    Sz_model_scatter = []
    Sz_ref_scatter = []
    s_m_max = S0 + alpha * math.pow((float(damelev) - Z0), beta)
    # TODO fix damelev with ZMAX when available
    s_r_max = ref_S0 + ref_alpha * math.pow((float(damelev) - ref_Z0), ref_beta)
    for h in z:
        s_m = S0 + alpha * math.pow((h - Z0), beta)
        Sz_model_scatter.append(s_m)

        s_r = ref_S0 + ref_alpha * math.pow((h - ref_Z0), ref_beta)
        Sz_ref_scatter.append(s_r*10000)

    Vs_model_scatter = []
    Vs_ref_scatter = []
    v_m_max = V0 + math.pow(((s_m_max-S0)/alpha), 1/beta) * (S0 + ((s_m_max-S0) / (beta + 1.)))
    v_r_max = ref_V0 + math.pow(((s_r_max-ref_S0)/ref_alpha), 1/ref_beta) * (ref_S0 + ((s_r_max-ref_S0) / (ref_beta + 1.)))
    Tx_model_scatter = []
    Tx_ref_scatter = []
    s_m_norm = []
    s_r_norm = []
    for a in s:
        v_m = V0 + math.pow(((a-S0)/alpha), 1/beta) * (S0 + ((a-S0) / (beta + 1.)))
        Vs_model_scatter.append(v_m)
        Tx_model_scatter.append(v_m/v_m_max)
        s_m_norm.append(a/s_m_max)

        #  v_r = ref_V0 + math.pow(((a*10000-ref_S0*10000)/ref_alpha), 1/ref_beta) * (ref_S0*10000 + ((a*10000-ref_S0*10000) / (ref_beta + 1.)))
        v_r = ref_V0 + math.pow(((a/10000-ref_S0)/ref_alpha), 1/ref_beta) * (ref_S0 + ((a/10000-ref_S0) / (ref_beta + 1.)))
        Vs_ref_scatter.append(v_r*10000)
        Tx_ref_scatter.append(v_r/v_r_max)
        s_r_norm.append(a/(s_r_max*10000))

    print(damelev)
    print(ref_Zmax)
    print(s_m_max)
    print(s_r_max*10000)
    print(v_m_max)
    print(v_r_max*10000)
    #  print(s)
    #  print(Vs_model_scatter)
    #  print(Vs_ref_scatter)

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
    Sz.set(xlabel='S(z)_reference (ha)',
           ylabel='S(z)_model (ha)')
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
    Sz_comp.set(xlabel='z (m)',
                ylabel='S(z) (ha)')
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
    Vs.set(xlabel='V(S)_reference (hm3) ',
           ylabel='V(S)_model (hm3)')
    # Trick to display in hm3
    #  Sz.xaxis.set_major_formatter(ticks_m3)
    #  Sz.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    Vs_comp = plt.subplot(gs[1])
    Vs_comp.plot(s, Vs_ref_scatter,
            'r-',
            label='V(S)_reference')
    Vs_comp.plot(s, Vs_model_scatter,
            'g-',
            label='V(S)_model')
    Vs_comp.axvline(x=s_r_max*10000,
                    ls=':', lw=2,
                    color='red',
                    label='Reference Max S')
    Vs_comp.axvline(x=s_m_max,
                    ls=':', lw=2,
                    color='green',
                    label='Model Max Dam S')
    Vs_comp.grid(b=True, which='major', linestyle='-')
    Vs_comp.grid(b=False, which='minor', linestyle='--')
    Vs_comp.set(xlabel='S(m2)',
           ylabel='V(S) (hm3)')
    # Trick to display in Ha
    #  Vs_comp.xaxis.set_major_formatter(ticks_m2)
    #  Vs_comp.yaxis.set_major_formatter(ticks_m3)
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
               label='Tx_model = f(Tx_reference)')
    Tx.grid(b=True, which='major', linestyle='-')
    Tx.grid(b=False, which='minor', linestyle='--')
    Tx.set(xlabel='Tx_reference (%) ',
           ylabel='Tx_model (%)')
    # Trick to display in hm3
    #  Sz.xaxis.set_major_formatter(ticks_m3)
    #  Sz.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    Tx_comp = plt.subplot(gt[1])
    Tx_comp.plot(s_r_norm, Tx_ref_scatter,
            'r-',
            label='Tx_reference')
    Tx_comp.plot(s_m_norm, Tx_model_scatter,
            'g-',
            label='Tx_model')
    Tx_comp.grid(b=True, which='major', linestyle='-')
    Tx_comp.grid(b=False, which='minor', linestyle='--')
    Tx_comp.set(xlabel='S/Smax (%)',
           ylabel='Tx (%)')
    # Trick to display in Ha
    #  Vs_comp.xaxis.set_major_formatter(ticks_m2)
    #  Vs_comp.yaxis.set_major_formatter(ticks_m3)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    plt.suptitle(damname + ": Tx",
                 fontsize=10)
    plt.savefig(os.path.splitext(args.outfile)[0]+"_Tx.png", dpi=300)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
