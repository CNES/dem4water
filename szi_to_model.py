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
import math
import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt
from statistics import median

def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i',
                        '--infile',
                        help="Input file")
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
        if ratio < 10:
            prev = sz
            stop_i = stop_i+1
        else:
            logging.debug("Dropping S_ZI after index "+ str(stop_i) +" with a delta ratio of "+str(ratio)+".")
            break_found=True
            break

    if break_found is True:
        Zi = Zi[stop_i:]
        S_Zi = S_Zi[stop_i:]

    Zi = Zi[::-1]
    S_Zi = S_Zi[::-1]
    logging.debug("Zi: ")
    logging.debug(Zi[:])
    logging.debug("S_Zi: ")
    logging.debug(S_Zi[:])
    logging.debug("Zi_max: "+ str(Zi[-1]) +" - S(Zi_max): "+ str(S_Zi[-1]))
    logging.debug("med(Zi):"+ str(median(Zi[1:])) +" - med(S_Zi): "+ str(median(S_Zi[1:])))
    logging.info("Z0: "+ str(Zi[0]) +" - S(Z0): "+ str(S_Zi[0]))

    #  # remove outliers / virtual surface overflow
    #  stop_i=0
    #  prev = S_Zi[1]
    #  for z, sz in zip(Zi, S_Zi):
    #      if prev != 0:
    #          ratio = sz/prev
    #      else:
    #          ratio = 1
    #      #  print(str(ratio))
    #      if ratio < 10:
    #          prev = sz
    #          stop_i = stop_i+1
    #      else:
    #          logging.debug("Dropping S_ZI after index "+ str(stop_i) +" with a delta ratio of "+str(ratio)+".")
    #          break
    #
    #  Zi = Zi[0:stop_i]
    #  S_Zi = S_Zi[0:stop_i]

    #  logging.debug("Zi: ")
    #  logging.debug(Zi[:])
    #  logging.debug("S_Zi: ")
    #  logging.debug(S_Zi[:])

    fig, ax = plt.subplots()
    ax.plot(Zi[1:], S_Zi[1:],
            marker='.',
            linestyle='dashed')
            #  markerfacecolor='blue')
    ax.set_ylim(bottom=-10.0)
    ax.set(xlabel='Virtual Water Surface Elevation (m)',
           ylabel='Virtual Water Surface (m2)')
    ax.plot(Zi[0], S_Zi[0],  'ro')
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=False, which='minor', linestyle='--')
    plt.minorticks_on()

    i=0
    best=-10000
    best_i=i
    best_alpha=0
    best_beta=0
    while ((i+10) < (len(Zi)-1)):
        #  p = np.polyfit(x[:-1], y[:-1], 1, rcond=None, full=False, w=None, cov=False)
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
                      + " and a local mae of: " +str(loc_mae))

        # Select MEA to be used:
        #  mae = glo_mae
        mae = loc_mae

        #  z = range(int(x[-1:])+1, int(x[0]))
        #  sz = []
        #  for h in z:
            #  s = y[-1:] + alpha * math.pow((h - x[-1:]), beta)
            #  sz.append(s)

        #  ax.plot(z, P(z), 'g+')
        #  ax.plot(z, sz, 'r--')

        #  if i == 54:
            #  break

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
    logging.info("S(Z) = "+format(S_Zi[0], '.2F')+" + "+format(best_alpha, '.3E')+ " * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'))
    plt.title("S(Z) = "+format(S_Zi[0], '.2F')+" + "+format(best_alpha, '.3E')+ " * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'))
    #  plt.title("S(Z) = "+str(S_Zi[0])+" + "+format(best_alpha, '.3E')+ " * ( Z - "+str(Zi[0])+" ) ^ "+ format(best_beta, '.3E'))

    z = range(int(Zi[0])+1, int(Zi[-1:]))
    sz = []
    for h in z:
        s = S_Zi[0] + best_alpha * math.pow((h - Zi[0]), best_beta)
        sz.append(s)

    #  ax.plot(z, P(z), 'g+')
    ax.plot(z, sz, 'r--')

    # Plot range / med(Zi); med(S_Zi)
    ax.plot(median(Zi[best_i:best_i+10]), median(S_Zi[best_i:best_i+10]), 'p', color='#ff7f0e')

    #  data_up = np.loadtxt(args.infile+"_up.dat")
    #  data_up = np.loadtxt(args.infile)
    #  x_up = data_up[:, 0]
    #  y_up = data_up[:, 1]
    #  plt.plot(x_up, y_up, 'go')

    #  plt.plot(x, y, 'ro')
    #  plt.plot(x[:-1], P(x[:-1]), 'go')
    #  plt.xlabel('Superficie (m²)')
    #  plt.ylabel('Cote du plan d\'eau (m)')


    #  ax.plot(x[:-1], P(x[:-1]), 'g+')
    #  plt.title(damname + ': S(Z_i)')
    plt.savefig(args.outfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
