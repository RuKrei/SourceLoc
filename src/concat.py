#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)


import os
import argparse
from os.path import join as opj
import glob
import mne
import matplotlib.pyplot as plt
import numpy as np
import logging


# Filter and resample
l_freq: float = 0.1    # lower pass-band edge
h_freq: float = 50.   # higher pass-band edge
fir_design: str = "firwin"
s_freq: int = 300



def define_subject(input_folder, args):
    subject = args.subject    
    if not subject:
        poss = [s for s in os.listdir(input_folder)]
        print(f"No subject specified, maybe you want to choose from those:\n {poss}")
        subject = input()
    if not subject.startswith("sub-"):
        ject = str(subject)
        subject = "sub-" + subject
    else:
        ject = subject.split("sub-")[-1]
    return subject, ject

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject", action="store", 
                        type=str, required=True,
                        help="Name of the Patient/ Subject to concatenate")
    parser.add_argument("--inputfolder", action="store", type=str, required=False, 
                        help="Specify a different data input folder")
    parser.add_argument("--fsalso", action="store", type=str, required=False, 
                        help="Use --fsonly true if you only want to do a freesurfer segmentation")      # do only freesurfer segmentation
    parser.add_argument("--openmp", action="store", type=str, required=False, 
                        help="Specify how many jobs/ processor cores to use")
    args = parser.parse_args()  
    
    # additional arguments
    if args.openmp:
        n_jobs = openmp = int(args.openmp)
    else:
        n_jobs = openmp = int(os.environ.get("OPENMP"))
        if n_jobs == None:
            n_jobs = openmp = int(1)

    if args.inputfolder:
        input_folder = args.inputfolder
    else:
        input_folder = os.environ.get("INPUT_FOLDER")

        
# define subject
    subject, ject = define_subject(input_folder, args)
        
# logging
    logfile = opj(input_folder, subject + "_concat.log")
    logging.basicConfig(filename=logfile, filemode="w",
                        format="\n%(levelname)s --> %(message)s")
    rootlog = logging.getLogger()
    rootlog.setLevel(logging.INFO)
    rootlog.info("Now running concatenation...")
    
# log parameters
    rootlog.info(f"*" * 20)
    rootlog.info("Parameters")
    rootlog.info(f"*" * 20)
    rootlog.info(f"Subject name = {ject}")
    rootlog.info(f"Input folder is set to: {input_folder}.")
    rootlog.info(f"Using {openmp} processor cores/ jobs.")
    rootlog.info(f"*" * 20)
    rootlog.info(f"Filtering files with margins {l_freq} and {h_freq}.")
    rootlog.info(f"Resampling frequency is set to {s_freq}")
    
