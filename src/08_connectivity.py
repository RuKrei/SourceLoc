#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sub", action="store", type=str, required=True)
args = parser.parse_args()

subj = args.sub


"""
To do:
generate a class (NxNConn) to compute NxN-Connectivity:
    inputs:
        fif
        metric (Coherence, PLI, GraphTheor.-Metrics)
            frequency and time-frequency metrics
        labels (--> subclass)
generate a class (OnexNConn) to compute 1xN-Connectivity:
    inputs:
        fif
        metric (Coherence, PLI, GraphTheor.-Metrics)
            frequency and time-frequency metrics
        source-seed/label
"""