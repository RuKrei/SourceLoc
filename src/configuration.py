#!/usr/bin/python

# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

subjects = ["BF28011991"]

#where am I?
beast =  True
candice = False

# BIDS inputs
if beast:
    bids_root = "/home/idrael/DATA/MEG/BIDS_clinic"   #beast
    data_root = "/home/idrael/DATA/MEG/clinic"    #beast

if candice:
    bids_root = "/Users/idrael/Playground/BIDS_clinic"    #candice
    data_root = "/Users/idrael/Playground/MEG/"           #candice
session = "resting"

# Processing
if beast:
    openmp = 35
    n_jobs = 35
if candice:
    openmp = 4
    n_jobs = 4

# data_file_preparation
concat_raws = True   # this only makes sense, if cHPI was on + Head position was transposed by Maxfilter
pick_meg = True    # Analysis in MEG
pick_eeg = False   # Analysis in EEG
# Filter
do_filter = True
l_freq = 1
h_freq = 70
fir_design = "firwin"
# Resample
do_resample = True   #downsampling jitters epoch-data. Don't overdo it + keep Nyquist in mind!
s_freq = 300


# freesurfer_and_BEM
do_anatomy = False
do_hippocampus_segmentation = True

# BEM_spacing = 
BEM_three_layer = [0.3, 0.006, 0.3]  # 3 layer BEM
BEM_single_shell = [0.3]
spacings = ["oct6", "oct5", "ico4", "ico5"]  # "oct5" = 1026, "ico4" = 2562, "oct6" = 4098, "ico5" = 10242


# Volume source space
volume_label = None    # standard freesurfer labels --> if set, then extra source spaces for the labels are constructed.
single_volume = True   # multiple values of volume label will be merged to a single source space, if true


# Frequency spectrum
freq_bands = dict(                #the frequency bands of interest for the analysis
                delta=(1, 4), 
                theta=(4, 7), 
                alpha=(8, 12), 
                beta=(13, 29), 
                gamma=(30, 80))
use_source_model_for_freq = "ico4"


# Source localization
use_source_model_for_sourceloc = "ico4"



"""
To do:
add a do_anatomy_only variable to do just the anatomy-stuff (or add a .sh file that does only that...)


"""