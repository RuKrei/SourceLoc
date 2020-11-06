#!/usr/bin/python

# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

subjects = ["PP05071984"]
do_anatomy = True


# BIDS inputs
bids_root = "/home/idrael/DATA/MEG/BIDS_clinic"
data_root = "/home/idrael/DATA/MEG_playground"
session = "resting"


# Processing
openmp = 35
n_jobs = 35
concat_raws = True   # this only makes sense, if cHPI was on + Head position was transposed by Maxfilter


# freesurfer_and_BEM
# BEM_spacing = 
BEM_three_layer = [0.3, 0.006, 0.3]  # 3 layer BEM
BEM_single_shell = [0.3]
spacings = ["oct6", "oct5", "ico4", "ico5"]  # "oct5" = 1026, "ico4" = 2562, "oct6" = 4098, "ico5" = 10242


# Volume source space
volume_label = None    # standard freesurfer labels --> if set, then extra source spaces for the labels are constructed.
single_volume = True   # multiple values of volume label will be merged to a single source space, if true