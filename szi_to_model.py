#!/usr/bin/env python3
'''
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
'''

""" szi_to_model.py
Prototype scrip allowing to derive a HSV model from a set of S(Z_i) values.
The first value of the set should be S_0 = S(Z_0) = 0 with Z_0 the altitude of dam bottom

The output
- the list of parameters [alpha;beta]
- a quality measurment of how well the model fit the S(Z_i) values
- additionnaly the plot of S(Z) and V(S)
"""

import os
import sys
import json
import math
import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from shapely.geometry import shape
from statistics import median


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i',
                        '--infile',
                        help="Input file")
    parser.add_argument('-d',
                        '--daminfo',
                        help="daminfo.json file")
    parser.add_argument('-z',
                        '--zmaxoffset',
                        type=int,
                        default=30,
                        help="Maximum elevation on top of dam elevation used for optimal model search")
    parser.add_argument('-o',
                        '--outfile',
                        help="Output file")
    parser.add_argument('--debug',
                        action='store_true',
                        help='Activate Debug Mode')


    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger('matplotlib').setLevel(logging.ERROR)

    args = parser.parse_args(arguments)

    print(args)

    logging_format = '%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s'
    if (args.debug is True):
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=logging_format)
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=logging_format)
    logging.info("Starting szi_to_model.py")

    # load GeoJSON file containing info
    with open(args.daminfo) as i:
        jsi = json.load(i)
    for feature in jsi['features']:
        if feature['properties']['name'] == 'Dam':
            logging.debug(feature)
            dam = shape(feature['geometry'])
            damname = feature['properties']['damname']
            damelev = feature['properties']['elev']


    data = np.loadtxt(args.infile)
    if (data.size <= 2):
        logging.error("Not enought S(Zi) data inside file "+args.infile)
        sys.exit("Error")

    Zi = data[:, 0]
    S_Zi = data[:, 1]

    # remove outliers / virtual surface overflow
    stop_i=0
    break_found=False
    prev = S_Zi[1]
    for z, sz in zip(Zi, S_Zi):
        if sz != 0:
            ratio = prev/sz
        else:
            ratio = 1
        #  print(str(ratio))
        #TODO: @parameters max ratio
        if ratio < 10:
            prev = sz
            stop_i = stop_i+1
        else:
            break_found=True
            break

    if break_found is True:
        logging.debug("Dropping S_ZI after index "+ str(stop_i) +" with a delta ratio of "+str(ratio)+".")
        Zi = Zi[stop_i:]
        S_Zi = S_Zi[stop_i:]
    else:
        logging.debug("No outliers detected, keeping all S_ZI data.")

    Zi = Zi[::-1]
    S_Zi = S_Zi[::-1]
    logging.debug("Zi: ")
    logging.debug(Zi[:])
    logging.debug("S_Zi: ")
    logging.debug(S_Zi[:])
    logging.debug("Zi_max: "+ str(Zi[-1]) +" - S(Zi_max): "+ str(S_Zi[-1]))
    logging.debug("med(Zi):"+ str(median(Zi[1:])) +" - med(S_Zi): "+ str(median(S_Zi[1:]))+ "(after outliers removal, Z0 and S(Z0) excluded)")
    logging.info("Z0: "+ str(Zi[0]) +" - S(Z0): "+ str(S_Zi[0]))

    i=0
    best=-10000
    best_i=i
    best_alpha=0
    best_beta=0
    l_mae=[]
    l_z=[]
    #TODO: @parameters winsize
    while ((i+10) < (len(Zi)-1)) and (median(Zi[i:i+10]) < args.zmaxoffset+float(damelev)):
        p = np.polyfit(Zi[i:i+10], S_Zi[i:i+10], 1, rcond=None, full=False, w=None, cov=False)
        P = np.poly1d(p)

        beta = p[0]*(median(Zi[i:i+10])-Zi[0])/(median(S_Zi[i:i+10])-S_Zi[0])
        alpha = p[0]*(math.pow(median(Zi[i:i+10])-Zi[0],1-beta))/beta

        mae_sum = 0
        for z, sz in zip(Zi, S_Zi):
            s = S_Zi[0] + alpha * math.pow((z - Zi[0]), beta)
            mae_sum += abs(sz - s)
        glo_mae = mae_sum / len(Zi)

        mae_sum = 0
        for z, sz in zip(Zi[i:i+10], S_Zi[i:i+10]):
            s = S_Zi[0] + alpha * math.pow((z - Zi[0]), beta)
            mae_sum += abs(sz - s)
        loc_mae = mae_sum / len(Zi[i:i+10])

        logging.debug('i: ' +str(i)
                      + " - Zrange [" +str(Zi[i])+ "; "+ str(Zi[i+10]) +"]"
                      + " - Slope: "+ str(p[0])
                      + " --> alpha= " +str(alpha)+ " - beta= " +str(beta)
                      + " with a glogal mae of: " +str(glo_mae)
                      + " m2 and a local mae of: " +str(loc_mae)
                      + " m2")

        # Select MEA to be used:
        #  mae = glo_mae
        mae = loc_mae
        l_z.append(median(Zi[i:i+10]))
        l_mae.append(mae)

        if (best == -10000) or (mae < best):
            best = mae
            best_i = i
            best_alpha = alpha
            best_beta = beta
            logging.debug("i: "+ str(i)
                          +" - alpha= " +str(alpha)
                          +" - beta= "  +str(beta)
                          +" - mae= " +str(mae))

        i = i+1


    logging.info('alpha= ' +str(best_alpha)+ " - beta= " +str(best_beta)
                 + " (Computed on the range ["
                 + str(Zi[best_i])+ "; "+ str(Zi[best_i+10]) +"])")
    logging.info(damname
                 +": S(Z) = "+format(S_Zi[0], '.2F')+" + "+format(best_alpha, '.3E')
                 +" * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'))
    z = range(int(Zi[0])+1, int(Zi[-1:]))
    sz = []
    for h in z:
        s = S_Zi[0] + best_alpha * math.pow((h - Zi[0]), best_beta)
        sz.append(s)

    # Moldel Plot
    fig, ax = plt.subplots()
    ax.plot(z, sz,
            'r--',
            label='S(z)')
    ax.plot(Zi[0], S_Zi[0],
            'ro')
    ax.scatter(Zi[1:], S_Zi[1:],
            marker='.',
            label='S(Zi)')
    ax.scatter(median(Zi[best_i:best_i+10]), median(S_Zi[best_i:best_i+10]),
               color='#ff7f0e',
               marker='p',
               label='Selected S(Zi)')
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=False, which='minor', linestyle='--')
    ax.set(xlabel='Virtual Water Surface Elevation (m)',
           ylabel='Virtual Water Surface (ha)')
    plt.title(damname
              +": S(Z) = "+format(S_Zi[0], '.2F')+" + "+format(best_alpha, '.3E')
              +" * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'),
              fontsize=10)
    ax.set_xlim(z[0]-10, z[-1]+10)
    # Trick to display in Ha
    ticks_m2 = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/10000.))
    ax.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')
    plt.savefig(args.outfile, dpi=300)

    # Superimpose local MAE to model
    #  ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    #  ax2.set_ylabel('Local MAE')  # we already handled the x-label with ax1
    #  ax2.set_yscale('log')
    #  ax2.yaxis.set_major_formatter(ticks_m2)
    #  ax2.plot(l_z, l_mae,
             #  marker='.',
             #  linestyle='dashed')
    #  ax2.plot(median(Zi[best_i:best_i+10]), best, 'p', color='#ff7f0e')
    #  plt.savefig(os.path.splitext(args.outfile)[0]+"_imposed.png", dpi=300)

    # Plot Local MAE
    maefig, maeax = plt.subplots()
    maeax.axvline(x=float(damelev),
                  ls='--', lw=1,
                  color='teal',
                  label='Dam elevation')
    maeax.plot(l_z, l_mae,
               marker='.',
               color='purple',
               linestyle='dashed',
               label='MAE')
    maeax.plot(median(Zi[best_i:best_i+10]), best,
               'p',
               color='#ff7f0e',
               label='min(MAE)')
    maeax.grid(b=True, which='major', linestyle='-')
    maeax.grid(b=False, which='minor', linestyle='--')
    maeax.set(xlabel='Virtual Water Surface Elevation (m)',
              ylabel='Local MAE (ha)')
    maeax.set_xlim(z[0]-10, z[-1]+10)
    maeax.set_yscale('log')
    # Trick to display in Ha
    maeax.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.title(damname+": Local Maximum Absolute Error")
    plt.legend(prop={'size': 6})
    if (args.debug is True):
        plt.savefig(os.path.splitext(args.outfile)[0]+"_mae.png", dpi=300)

    # Combined Local MAE / model plot
    fig = plt.figure(dpi=300)
    fig.subplots_adjust(hspace=0)
    gs = fig.add_gridspec( 2, 1, height_ratios=[4, 1])

    maeax0 = plt.subplot(gs[0])
    maeax0.axvline(x=float(damelev),
                   ls=':', lw=2,
                   color='teal',
                   label='Dam Elevation')
    maeax0.plot(z, sz,
                'r--',
                label='S(z)')
    maeax0.plot(Zi[0], S_Zi[0],
                'ro')
    maeax0.scatter(Zi, S_Zi,
                   marker='.',
                   label='S(Z_i)')
    maeax0.plot(Zi[best_i:best_i+10], S_Zi[best_i:best_i+10],
                marker='o',
                ls=':', lw=1,
                color='yellow',
                markerfacecolor='blue',
                label='Selected S(Z_i)')
    maeax0.plot(median(Zi[best_i:best_i+10]), median(S_Zi[best_i:best_i+10]),
                color='#ff7f0e',
                marker='p',
                markerfacecolor='#ff7f0e')
    maeax0.grid(b=True, which='major', linestyle='-')
    maeax0.grid(b=False, which='minor', linestyle='--')
    maeax0.set(xlabel='Virtual Water Surface Elevation (m)',
               ylabel='Virtual Water Surface (ha)')
    maeax0.label_outer()
    maeax0.set_xlim(z[0]-10, z[-1]+10)
    # Trick to display in Ha
    maeax0.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='upper left')

    maeax1 = plt.subplot(gs[1])
    maeax1.axvline(x=float(damelev),
                   ls=':', lw=2,
                   color='teal',
                   label='Dam Elevation')
    maeax1.plot(l_z, l_mae,
                color='purple',
                marker='.',
                linestyle='dashed',
                label='MAE')
    maeax1.plot(median(Zi[best_i:best_i+10]), best,
                'p',
                color='#ff7f0e',
                label='min(MAE)')
    maeax1.grid(b=True, which='major', linestyle='-')
    maeax1.grid(b=True, which='minor', linestyle='--')
    maeax1.set(xlabel='Virtual Water Surface Elevation (m)',
               ylabel='Local MAE (ha)')
    maeax1.set_xlim(z[0]-10, z[-1]+10)
    maeax1.set_yscale('log')
    # Trick to display in Ha
    maeax1.yaxis.set_major_formatter(ticks_m2)
    plt.minorticks_on()
    plt.legend(prop={'size': 6}, loc='lower left')

    plt.suptitle(damname + ": S(Z) = "+format(S_Zi[0], '.2F')
                 +" + "+format(best_alpha, '.3E')
                 +" * ( Z - "+format(Zi[0], '.2F')
                 +" ) ^ "+ format(best_beta, '.3E'),
                 fontsize=10)
    plt.savefig(os.path.splitext(args.outfile)[0]+"_combo.png", dpi=300)

    # V(S)
    S_dam = S_Zi[0] + best_alpha * math.pow((float(damelev) - Zi[0]), best_beta)
    s = range(0, math.ceil(S_dam)+10000)
    v0 = 0.
    vs = []
    for a in s:
        v = v0 + math.pow(((a-S_Zi[0])/best_alpha), 1/best_beta) * (S_Zi[0] + ((a-S_Zi[0]) / (best_beta + 1.)))
        vs.append(v)

    #  # V(S) Moldel Plot
    #  figv, axv = plt.subplots()
    #  axv.axvline(x=S_dam,
    #              ls=':', lw=2,
    #              color='teal',
    #              label='Max Dam Surface')
    #  axv.plot(s, vs,
    #          'r-',
    #          label='V(S)')
    #  axv.grid(b=True, which='major', linestyle='-')
    #  axv.grid(b=False, which='minor', linestyle='--')
    #  axv.set(xlabel='Estimated Water Surface (ha)',
    #          ylabel='Modelized Water Volume (hm3)')
    #  plt.title(damname + " : V=f(S) ",
    #            #  +": S(Z) = "+format(S_Zi[0], '.2F')+" + "+format(best_alpha, '.3E')
    #            #  +" * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'),
    #            fontsize=10)
    #  # Trick to display in Ha
    #  axv.xaxis.set_major_formatter(ticks_m2)
    #  ticks_m3 = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000000.))
    #  axv.yaxis.set_major_formatter(ticks_m3)
    #  plt.minorticks_on()
    #  plt.legend(prop={'size': 6}, loc='upper left')
    #  plt.savefig(os.path.splitext(args.outfile)[0]+"_VS.png", dpi=300)



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
