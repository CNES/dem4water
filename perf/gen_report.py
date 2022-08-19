#!/usr/bin/env python3
"""
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CNES
:copyright: 2022 CNES. All rights reserved.
:license: see LICENSE file
:created: 2022
"""

import logging
import subprocess
import sys


def get_current_git_rev():
    return (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("ascii")
        .strip()
    )


def main(arguments):
    """gen_report.py
    Entry point to generate performance report of the current revision
    """

    # Setup Logger
    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=logging_format)

    # Expose revision
    logging.info("Performance Estimation of revision " + get_current_git_rev())

    # Collect baseline sites from parameters


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
