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
    parser.add_argument('--loglevel',
                        help="Log Level",
                        default="INFO")

    args = parser.parse_args(arguments)

    print(args)

    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.loglevel)
    logging_format = '%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s'
    logging.basicConfig(stream=sys.stdout, level=numeric_level, format=logging_format)
    logging.info("Starting szi_to_model.py")

    data = np.loadtxt(args.infile)
    x = data[:, 0]
    y = data[:, 1]

    print("x: ")
    print(x[:-1])

    logging.info("med(zi):"+ str(median(x[:-1])) +" - med(Szi): "+ str(median(y[:-1])))

    plt.plot(x, y, 'ro')
    plt.savefig(args.outfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
