#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from genericpath import isdir
from os.path import isfile
from nipype.interfaces.freesurfer import ReconAll
import os
from configuration import (openmp, n_jobs, do_anatomy, bids_root, data_root, do_hippocampus_segmentation,
                            BEM_single_shell, BEM_three_layer, spacings, volume_label, single_volume, 
                            derivatives_root, subjects_dir)
import glob
import mne
from utils.utils import FileNameRetriever, recursive_overwrite
import shutil
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sub", action="store", type=str, required=True)
args = parser.parse_args()

subj = args.sub

fnr = FileNameRetriever(derivatives_root)

def run_shell_command(command):
    subprocess.run(command, shell=True, capture_output=True, check=True)


# 1. freesurfer
if do_anatomy == True:
    reconall = ReconAll()
    if not subj.startswith("sub-"):
        subj = "sub-" + subj
    # check if freesurfer segmentation was already performed
    this_subjects_dir = fnr.get_filename(subj, "subjects_dir")
    freesurfered = os.path.join(subjects_dir, subj)
    if not isdir(freesurfered):
        anafolder = os.path.join(bids_root, subj, "ses-resting", "anat")
        nii_file = glob.glob(anafolder + "/*.nii*")[0]
        subjects_dir = subjects_dir
        reconall.inputs.subject_id = subj
        reconall.inputs.T1_files = nii_file
        reconall.inputs.directive = 'all'
        reconall.inputs.subjects_dir = subjects_dir
        reconall.inputs.openmp = openmp
        reconall.inputs.flags = "-3T"
        reconall.run()
    else:
        print(f"A freesurfer segmentation of subject {subj} already exists in {subjects_dir}")
        freesurfered = os.path.join(subjects_dir, subj)
    if do_hippocampus_segmentation:
        hippofile = os.path.join(freesurfered, "mri", "lh.hippoSfVolumes-T1.v21.txt")
        if not isfile(hippofile):
            try:
                print(f"Now running hippocampal segmentation for subject: {subj}\nThis might take some time")
                hipposeg = "segmentHA_T1.sh " + subj
                run_shell_command(hipposeg)
            except Exception as e:
                print(f"Something went wrong with the hippocampal segmentation! --> {e}")
        else:
            print(f"Omitting hippocampal segmentation for subject {subj}, as it already exists")
       
    local_copy = os.path.join(this_subjects_dir, subj)
    if not isdir(local_copy):
        print(f"Copying freesufer folder of {subj} to {this_subjects_dir}")
        recursive_overwrite(freesurfered, local_copy)



# from here on subjects_dir refers to the local subjects_dir in .derivatives/subj/freesurfer - folder

if not subj.startswith("sub-"):
    subj = "sub-" + subj
this_subjects_dir = fnr.get_filename(subj, "subjects_dir")


# Head-Model
try:
    if not os.path.isfile(this_subjects_dir + "/" + subj + "/bem/" + subj + "-head.fif"):
        mne.bem.make_watershed_bem(subject=subj, subjects_dir=this_subjects_dir, overwrite=True)
except Exception as e:
    print("#"*30)
    print("#"*30)
    print("#"*30)
    print(f"Failed to make watershed BEM for {subj} --> {e}")


# Cortex Source space
for spacing in spacings:
    srcfilename = fnr.get_filename(subj, spacing)
    if not os.path.isfile(srcfilename):
        try:
            src = mne.setup_source_space(subj, spacing = spacing, 
                                            subjects_dir = this_subjects_dir, 
                                            n_jobs=n_jobs, 
                                            verbose=True)
            mne.write_source_spaces(srcfilename, src, overwrite=True, verbose=True)
        except Exception as e:
            print("#"*30)
            print("#"*30)
            print("#"*30)
            print(f"Failed to setup source space with spacing {spacing} for {subj} --> {e}")


# Volume source space
srcfilename = fnr.get_filename(subj, "vol-src")
if not os.path.isfile(srcfilename):
    try:    
        src_vol = mne.setup_volume_source_space(subj, pos=3.0, 
                                        subjects_dir = this_subjects_dir, 
                                        volume_label=volume_label,
                                        single_volume=single_volume,
                                        verbose=True)
        mne.write_source_spaces(srcfilename, src_vol, overwrite=True, verbose=True)
    except Exception as e:
        print("#"*30)
        print("#"*30)
        print("#"*30)
        print(f"Failed to setup Volume source space for {subj} --> {e}")


# BEM Solutions - single shell
bem_save_name = fnr.get_filename(subj, "single-shell-model")
if not os.path.isfile(bem_save_name):
    try:
        bem = mne.make_bem_model(subj, ico=4, 
                        conductivity=BEM_single_shell,   
                        subjects_dir=this_subjects_dir, verbose=True)
        mne.write_bem_surfaces(bem_save_name, bem, overwrite=True) #, overwrite=True)
    except Exception as e:
        print("#"*30)
        print("#"*30)
        print("#"*30)
        print(f"Failed to setup single shell BEM model for {subj} --> {e}")
bem_sol_filename = fnr.get_filename(subj, "single-shell-BEM-sol")
if not os.path.isfile(bem_sol_filename):
    try:
        bem = mne.read_bem_surfaces(bem_save_name)    
        bem_sol = mne.make_bem_solution(bem)
        mne.write_bem_solution(bem_sol_filename, bem_sol, overwrite=True)
    except Exception as e:
        print("#"*30)
        print("#"*30)
        print("#"*30)
        print(f"Failed to calculate BEM solution (single-shell) for {subj} --> {e}")


# BEM Solutions - 3-layer-BEM
bem_save_name = fnr.get_filename(subj, "3-layer-BEM-model")
if not os.path.isfile(bem_save_name):
    try:
        bem = mne.make_bem_model(subj, ico=4, 
                        conductivity=BEM_three_layer,   
                        subjects_dir=this_subjects_dir, verbose=True)
        mne.write_bem_surfaces(bem_save_name, bem, overwrite=True) #, overwrite=True)
    except Exception as e:
        print("#"*30)
        print("#"*30)
        print("#"*30)
        print(f"Failed to calculate 3-layer BEM model for {subj} --> {e}")
bem_sol_filename = fnr.get_filename(subj, "3-layer-BEM-sol")
if not os.path.isfile(bem_sol_filename):
    try:
        bem = mne.read_bem_surfaces(bem_save_name)    
        bem_sol = mne.make_bem_solution(bem)
        mne.write_bem_solution(bem_sol_filename, bem_sol, overwrite=True)
    except Exception as e:
        print("#"*30)
        print("#"*30)
        print("#"*30)
        print(f"Failed to calculate 3-layer BEM solution for {subj} --> {e}")
        print("This is bad, please look into the freesurfer segmentation...")
        print("Alternatively, you might be able to run the analysis with a single-shell-head-model (look into the configuration file")

