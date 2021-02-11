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

    args = parser.parse_args(arguments)

    print(args)

    logging_format = '%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s'
    if (args.debug is True):
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=logging_format)
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=logging_format)
    logging.info("Starting szi_to_model.py")

    data = np.loadtxt(args.infile)
    x = data[:, 0]
    y = data[:, 1]
    print("x: ")
    print(x[:-1])
    logging.info("med(zi):"+ str(median(x[:-1])) +" - med(Szi): "+ str(median(y[:-1])))
    logging.info("z0:"+ str(x[-1:]) +" - S0: "+ str(y[-1:]))

    i=0
    while ((i+10) < (len(x)-1)):
        #  p = np.polyfit(x[:-1], y[:-1], 1, rcond=None, full=False, w=None, cov=False)
        p = np.polyfit(x[i:i+10], y[i:i+10], 1, rcond=None, full=False, w=None, cov=False)
        logging.info('i: ' +str(i)+ " - range [" +str(x[i])+ "; "+ str(x[i+10]) +"]")
        print(p)
        beta = p[0]*(median(x[i:i+10])-x[-1:])/(median(y[i:i+10])-y[-1:])
        alpha = p[0]*(math.pow(median(x[i:i+10])-x[-1:],1-beta))/beta
        logging.info('alpha= ' +str(alpha)+ " - beta= " +str(beta))

        if i == 54:
            break

        i = i+1
        P = np.poly1d(p)

    z = range(int(x[-1:])+1, int(x[0]))
    sz = []
    for h in z:
        s = y[-1:] + alpha * math.pow((h - x[-1:]), beta)
        sz.append(s)

    #  data_up = np.loadtxt(args.infile+"_up.dat")
    #  data_up = np.loadtxt(args.infile)
    #  x_up = data_up[:, 0]
    #  y_up = data_up[:, 1]
    #  plt.plot(x_up, y_up, 'go')


    #  plt.plot(x, y, 'ro')
    #  plt.plot(x[:-1], P(x[:-1]), 'go')
    #  plt.xlabel('Superficie (m²)')
    #  plt.ylabel('Cote du plan d\'eau (m)')


    fig, ax = plt.subplots()
    #  ax.plot(x[:-1], y[:-1],
    ax.plot(x[:-1], y[:-1],
            color='#ff7f0e',
            marker='o',
            linestyle='dashed',
            markerfacecolor='blue')
    #  ax.plot(x[:-1], P(x[:-1]), 'g+')
    ax.plot(z, P(z), 'g+')
    ax.plot(z, sz, 'r--')
    ax.set_ylim(bottom=-10.0)
    ax.set(xlabel='Virtual Water Surface Elevation (m)',
           ylabel='Virtual Water Surface (m2)')
    ax.plot(x[-1], y[-1],  'ro')
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=False, which='minor', linestyle='--')
    plt.minorticks_on()
    #  plt.title(damname + ': S(Z_i)')
    plt.savefig(args.outfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
