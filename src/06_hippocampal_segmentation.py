# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from nipype.interfaces.freesurfer import ReconAll, FSCommand
import os
from configuration import (subjects, openmp, n_jobs, do_anatomy, bids_root, data_root, 
                        do_hippocampus_segmentation)
import glob
from utils.utils import FileNameRetriever

                            
fnr = FileNameRetriever(bids_root)

if do_hippocampus_segmentation:
    for subj in subjects:  
        anafolder = os.path.join(data_root, subj, "data", "anat", subj)
        nii_file = glob.glob(anafolder + "/*.nii*")
        subjects_dir = fnr.get_filename(subj, file="subjects_dir")
        reconall_subfields = ReconAll()
        reconall_subfields.inputs.subject_id = subj
        reconall_subfields.inputs.directive = 'all'
        reconall_subfields.inputs.subjects_dir = subjects_dir
        reconall_subfields.inputs.T1_files = nii_file
        reconall_subfields.inputs.hippocampal_subfields_T1 = True
        reconall_subfields.inputs.openmp = openmp
        reconall_subfields.cmdline




"""
freesurfer - hippocampal segmentation, cite:
Iglesias, J.E., Augustinack, J.C., Nguyen, K., Player, C.M., Player, A., Wright,
M., Roy, N., Frosch, M.P., McKee, A.C., Wald, L.L., Fischl, B., and Van Leemput, K.,
A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution
MRI: Application to adaptive segmentation of in vivo MRI.  Neuroimage 115, 2015, 117-137.
http://dx.doi.org/10.1016/j.neuroimage.2015.04.042

"""



"""
To do:
use freesurfers segmentHA_T1.sh $subject script
visualize single subjects segmentation and statistics:
    right to left difference?
    hippocampi too small relative to brain volume?
visualize single subject and statistics vs. grand average of healthy subjects in same age bin

"""