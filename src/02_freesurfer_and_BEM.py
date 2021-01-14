#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from genericpath import isdir
from os.path import isfile
from nipype.interfaces.freesurfer import ReconAll
import os
from configuration import (subjects, openmp, n_jobs, do_anatomy, bids_root, data_root, do_hippocampus_segmentation,
                            BEM_single_shell, BEM_three_layer, spacings, volume_label, single_volume, derivatives_root,
                            subjects_dir)
import glob
import mne
from utils.utils import FileNameRetriever
import shutil
import subprocess


fnr = FileNameRetriever(derivatives_root)

def recursive_overwrite(src, dest, ignore=None):
    if os.path.isdir(src):
        if not os.path.isdir(dest):
            os.makedirs(dest)
        files = os.listdir(src)
        if ignore is not None:
            ignored = ignore(src, files)
        else:
            ignored = set()
        for f in files:
            if f not in ignored:
                recursive_overwrite(os.path.join(src, f), 
                                    os.path.join(dest, f), 
                                    ignore)
    else:
        shutil.copyfile(src, dest)

def run_shell_command(command):
    subprocess.run(command, shell=True, capture_output=True, check=True)


# 1. freesurfer
if do_anatomy == True:
    reconall = ReconAll()
    for subj in subjects:
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

        if do_hippocampus_segmentation:
            hippofile = os.path.join(freesurfered, "mri", "lh.hippoSfVolumes-T1.v21.txt")
            if not isfile(hippofile):
                print(f"Now running hippocampal segmentation for subject: {subj}\nThis might take some time")
                hipposeg = "segmentHA_T1.sh " + subj
                run_shell_command(hipposeg)
            else:
                print(f"Omitting hippocampal segmentation for subject {subj}, as it already exists")
           
        local_copy = os.path.join(this_subjects_dir, subj)
        if not isdir(local_copy):
            print(f"Copying freesufer folder of {subj} to {this_subjects_dir}")
            recursive_overwrite(freesurfered, local_copy)


# from here on subjects_dir refers to the local subjects_dir in .derivatives/subj/freesurfer - folder

for subj in subjects:
    if not subj.startswith("sub-"):
        subj = "sub-" + subj

# Head-Model
    subjects_dir = fnr.get_filename(subj, file="subjects_dir")
    if not os.path.isfile(subjects_dir + "/" + subj + "/bem/" + subj + "-head.fif"):
        mne.bem.make_watershed_bem(subject=subj, subjects_dir=subjects_dir, overwrite=True)


# Cortex Source space
    for spacing in spacings:
        srcfilename = fnr.get_filename(subj, spacing)
        if not os.path.isfile(srcfilename):
            src = mne.setup_source_space(subj, spacing = spacing, 
                                            subjects_dir = subjects_dir, 
                                            n_jobs=n_jobs, 
                                            verbose=True)
            mne.write_source_spaces(srcfilename, src, overwrite=True, verbose=True)


# Volume source space
    srcfilename = fnr.get_filename(subj, "vol-src")
    if not os.path.isfile(srcfilename):
        src_vol = mne.setup_volume_source_space(subj, pos=3.0, 
                                        subjects_dir = subjects_dir, 
                                        volume_label=volume_label,
                                        single_volume=single_volume,
                                        verbose=True)
        mne.write_source_spaces(srcfilename, src_vol, overwrite=True, verbose=True)


# BEM Solutions - single shell
    bem_save_name = fnr.get_filename(subj, "single-shell-model")
    bem = mne.make_bem_model(subj, ico=4, 
                    conductivity=BEM_single_shell,   
                    subjects_dir=subjects_dir, verbose=True)
    mne.write_bem_surfaces(bem_save_name, bem, overwrite=True) #, overwrite=True)

    bem_sol_filename = fnr.get_filename(subj, "single-shell-BEM-sol")
    if not os.path.isfile(bem_sol_filename):
        bem = mne.read_bem_surfaces(bem_save_name)    
        bem_sol = mne.make_bem_solution(bem)
        mne.write_bem_solution(bem_sol_filename, bem_sol, overwrite=True)


# BEM Solutions - 3-layer-BEM
    bem_save_name = fnr.get_filename(subj, "3-layer-BEM-model")
    bem = mne.make_bem_model(subj, ico=4, 
                    conductivity=BEM_three_layer,   
                    subjects_dir=subjects_dir, verbose=True)
    mne.write_bem_surfaces(bem_save_name, bem, overwrite=True) #, overwrite=True)

    bem_sol_filename = fnr.get_filename(subj, "3-layer-BEM-sol")
    if not os.path.isfile(bem_sol_filename):
        bem = mne.read_bem_surfaces(bem_save_name)    
        bem_sol = mne.make_bem_solution(bem)
        mne.write_bem_solution(bem_sol_filename, bem_sol, overwrite=True)

