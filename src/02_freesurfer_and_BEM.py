# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from nipype.interfaces.freesurfer import ReconAll
import os
from configuration import (subjects, openmp, n_jobs, do_anatomy, bids_root, data_root,
                            BEM_single_shell, BEM_three_layer, spacings, volume_label, single_volume)
import glob
import mne
from utils.utils import FileNameRetriever


# freesurfer
if do_anatomy == True:
    reconall = ReconAll()
    reconall_subfields = ReconAll()
    for subj in subjects:
        anafolder = os.path.join(data_root, subj, "data", "anat", subj)
        nii_file = glob.glob(anafolder + "/*.nii*")
        subjects_dir = os.path.join(bids_root, "derivatives", "sub-" + subj, "freesurfer")
        reconall.inputs.subject_id = subj
        reconall.inputs.T1_files = nii_file
        reconall.inputs.directive = 'all'
        #reconall.inputs.hippocampal_subfields_T1 = True
        reconall.inputs.subjects_dir = subjects_dir
        reconall.inputs.openmp = openmp
        reconall.run()

        
        


fnr = FileNameRetriever(bids_root)


# Head-Model
for subj in subjects:
    subjects_dir = fnr.get_filename(subj, file="subjects_dir")
    if not os.path.isfile(subjects_dir + '/' + subj + '/bem/' + subj + '-head.fif'):
        mne.bem.make_watershed_bem(subject=subj, subjects_dir=subjects_dir, overwrite=True)


# Cortex Source space
for subj in subjects:
    for spacing in spacings:
        srcfilename = fnr.get_filename(subj, spacing)
        subjects_dir = fnr.get_filename(subj, file="subjects_dir")
        if not os.path.isfile(srcfilename):
            src = mne.setup_source_space(subj, spacing = spacing, 
                                            subjects_dir = subjects_dir, 
                                            n_jobs=n_jobs, 
                                            verbose=True)
            mne.write_source_spaces(srcfilename, src, overwrite=True, verbose=True)


# Volume source space
for subj in subjects:
    srcfilename = fnr.get_filename(subj, "vol-src")
    subjects_dir = fnr.get_filename(subj, file="subjects_dir")
    if not os.path.isfile(srcfilename):
        src_vol = mne.setup_volume_source_space(subj, pos=3.0, 
                                        subjects_dir = subjects_dir, 
                                        volume_label=volume_label,
                                        single_volume=single_volume,
                                        verbose=True)
        mne.write_source_spaces(srcfilename, src_vol, overwrite=True, verbose=True)


# BEM Solutions - single shell
for subj in subjects:
    subjects_dir = fnr.get_filename(subj, file="subjects_dir")
    bem_save_name = fnr.get_filename(subj, "single-shell-model")
    bem = mne.make_bem_model(subj, ico=4, 
                    conductivity=BEM_single_shell,   
                    subjects_dir=subjects_dir, verbose=True)
    mne.write_bem_surfaces(bem_save_name, bem) #, overwrite=True)

    bem_sol_filename = fnr.get_filename(subj, "single-shell-BEM-sol")
    if not os.path.isfile(bem_sol_filename):
        bem = mne.read_bem_surfaces(bem_save_name)    
        bem_sol = mne.make_bem_solution(bem)
        mne.write_bem_solution(bem_sol_filename, bem_sol)


# BEM Solutions - 3-layer-BEM
for subj in subjects:
    subjects_dir = fnr.get_filename(subj, file="subjects_dir")
    bem_save_name = fnr.get_filename(subj, "3-layer-BEM-model")
    bem = mne.make_bem_model(subj, ico=4, 
                    conductivity=BEM_three_layer,   
                    subjects_dir=subjects_dir, verbose=True)
    mne.write_bem_surfaces(bem_save_name, bem) #, overwrite=True)

    bem_sol_filename = fnr.get_filename(subj, "3-layer-BEM-sol")
    if not os.path.isfile(bem_sol_filename):
        bem = mne.read_bem_surfaces(bem_save_name)    
        bem_sol = mne.make_bem_solution(bem)
        mne.write_bem_solution(bem_sol_filename, bem_sol)

