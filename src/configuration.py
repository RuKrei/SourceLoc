#!/usr/bin/python

# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)


import os

subjects = ["GM30091987", "ARM23082007"]                 # ["AD14071990", "PW17081978"]   # "HA11041987",  #["KF28091994"]  # ["FT05042011"] 
                            #["LL23052000"]  # ["BF28011991"] #  ["MKE03101965"]  
                            # ["johndoe"] # ["MA17011989"]   # ["PP05071984"] #  # "*" means every subject in data_root

subjects_dir = os.environ.get("SUBJECTS_DIR")

#where am I?
beast =  True
candice = False
h_beast = False

# BIDS inputs
if beast:
    bids_root = "/home/idrael/DATA/MEG/SourceLocTest/BIDSTestData"      #beast
    data_root = "/home/idrael/DATA/MEG/SourceLocTest/input_folder"      #beast
    input_folder = "/home/idrael/DATA/MEG/new_patients/"               #change, if data should come from another folder

if candice:
    bids_root = "/Users/idrael/Playground/SourceLocTest/BIDSTestData"           #candice
    data_root = "/Users/idrael/Playground/SourceLocTest/input_folder"           #candice

if h_beast:
    bids_root = "/media/idrael/Data/SourceLocTest/BIDSTestData"         #h_beast
    data_root = "/media/idrael/Data/SourceLocTest/input_folder"         #h_beast
    subjects_dir = "/media/idrael/Data/Playground/anat"                 #h_beast

session = "resting"
derivatives_root = os.path.join(bids_root, "derivatives")
extras_dir = os.path.join(data_root, "extras")

if subjects == ["*"]:
    subjects = os.listdir(data_root)

# Processing
if beast:
    openmp = 35
    n_jobs = 35
if candice or h_beast:
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
# Artifacts
n_grad = 1          # number of Vectors to apply for gradiometers
n_mag = 1           # number of Vectors to apply for magnetometers
n_eeg = 1           # number of Vectors to apply for eeg
ecg_channel=None
do_ecg_correction_ssp = True
do_ecg_correction_ica = False
do_ecg_correction_regression = False

eog_channel=None
do_eog_correction_ssp = False
do_eog_correction_ica = False
do_eog_correction_regression = False


# freesurfer_and_BEM
do_anatomy = True
do_hippocampus_segmentation = True


# Coregistration
use_single_transfile = True     # the pipeline assumes, that .fifs have been maxfiltered with head translation (--trans option)
                                # that way one trans-file should suffice

# BEM_spacing = 
BEM_three_layer = [0.3, 0.006, 0.3]  # 3 layer BEM
BEM_single_shell = [0.3]
spacings = ["oct6", "oct5", "ico4", "ico5"]  # "oct5" = 1026, "ico4" = 2562, "oct6" = 4098, "ico5" = 10242


# Volume source space
volume_label = None    # standard freesurfer labels --> if set, then extra source spaces for the labels are constructed (to do).
single_volume = True   # multiple values of volume label will be merged to a single source space, if true


# Frequency spectrum
freq_bands = dict(                #the frequency bands of interest for the analysis
                delta=(1, 4), 
                theta=(4, 7), 
                alpha=(8, 12), 
                beta=(13, 29), 
                gamma=(30, 80))
use_source_model_for_freq = "ico4"


# Source localization                                   # currently only works with use_single_transfile = True 
#                                                         (see configuration of coregistration)
do_source_loc = True                                    # perform source localization fro events

use_source_model_for_sourceloc = "ico4"                 # Options = "oct6", "oct5", "ico4", "ico5"
do_volume_source_loc = False
use_fwd_model_for_sourceloc = "single-shell-BEM-sol" #"3-layer-BEM-sol"         # or:   "single-shell-BEM-sol"  
source_loc_methods = ["dSPM", "eLORETA"]                # "MNE"
vol_source_loc_methods = "eLORETA"
signal_to_noise_ratio = 3.
minimum_norm_ori = None                                 # None --> norm of loose/free orientations; 
                                                        # "normal" --> only 90° to surface
inv_loose_option = 1.                                   # 1. == free; 0.2 == constrained.
peaks_tmin = -.025
peaks_tmax = 0.
peaks_mode = "abs"                                      # How to deal with the sign of the data:  "pos", "neg", or "abs"
peaks_nr_of_points = 5                                   

dip_times = {   'min20ms':  (-0.025,-0.020),            # Time points in Miliseconds for Equivalent current dipole fit
                'min15ms':  (-0.019,-0.015),
                'min10ms':  (-0.014,-0.010),
                'min5ms':  (-0.009,-0.005),
                'peak':     (-0.004,0.000)}
