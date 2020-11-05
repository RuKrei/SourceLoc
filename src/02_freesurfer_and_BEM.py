# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from nipype.interfaces.freesurfer import ReconAll
import os
from configuration import (subjects, openmp, n_jobs, do_anatomy, bids_root, data_root,
                            BEM_single_shell, BEM_three_layer, spacings)
import glob
import mne


# freesurfer
reconall = ReconAll()
for subj in subjects:
    anafolder = os.path.join(data_root, subj, "data", "anat", subj)
    nii_file = glob.glob(anafolder + "/*.nii*")
    subjects_dir = os.path.join(bids_root, "derivatives", "sub-" + subj, "freesurfer")
    reconall.inputs.subject_id = subj
    reconall.inputs.T1_files = nii_file
    reconall.inputs.directive = 'all'
    reconall.inputs.subjects_dir = subjects_dir
    reconall.inputs.openmp = openmp
    reconall.run()


# Cortex Source space
for subj in subjects:
    for spacing in spacings:
        fbase = os.path.join(bids_root, "derivatives", "sub-" + subj)
        fsrc = os.path.join(fbase, "source_model")
        srcfilename = (subj + '_' + spacing + '_src.fif')
        srcfilename = os.path.join(fsrc, srcfilename)
        if not os.path.isfile(srcfilename):
                src = mne.setup_source_space(subj, spacing = spacing, 
                                            subjects_dir = subjects_dir, 
                                            n_jobs=n_jobs, 
                                            verbose=True)
                mne.write_source_spaces(srcfilename, src, overwrite=True, verbose=True)


# BEM-Model
#for subj in subjects:
#    subjects_dir = os.path.join(bids_root, "derivatives", "sub-" + subj, "freesurfer")
#    if not os.path.isfile(subjects_dir + '/' + subj + '/bem/' + subj + '-head.fif'):
#        mne.bem.make_watershed_bem(subject=subj, subjects_dir=subjects_dir, overwrite=True)
    






# BEM Solutions    
    
#bem_f_name = 



#    if not os.path.isfile(bem_f_name):
#        bem = mne.make_bem_model(subj, ico=bem_spacing, 
#                        conductivity=conductivity,   
#                        subjects_dir=subjects_dir, verbose=True)
#        bem_save_name = os.path.join(bem_f_name)
#        mne.write_bem_surfaces(bem_save_name, bem)
#
#    bem_sol_filename = (subj + '_bem.fif')
#    bem_sol_filename = os.path.join(for_report_dir, bem_sol_filename)
#    if not os.path.isfile(bem_sol_filename):
#        if not bem:
#            bem = mne.read_bem_surfaces(bem_file)    
#        bem_sol = mne.make_bem_solution(bem)
#        mne.write_bem_solution(bem_sol_filename, bem_sol)
#    else:
#        bem_sol = mne.bem.read_bem_solution(bem_sol_filename)