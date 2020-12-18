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


    data_up = np.loadtxt(args.infile+"_up.dat")
    x_up = data_up[:, 0]
    y_up = data_up[:, 1]


    plt.plot(x, y, 'ro')
    plt.plot(x_up, y_up, 'go')
    plt.xlabel('Superficie (m²)')
    plt.ylabel('Cote du plan d\'eau (m)')
    plt.savefig(args.outfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
