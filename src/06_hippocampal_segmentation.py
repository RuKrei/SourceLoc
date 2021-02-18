# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from nipype.interfaces.freesurfer import ReconAll, FSCommand
import os
from configuration import (subjects, openmp, n_jobs, do_anatomy, bids_root, data_root, 
                        do_hippocampus_segmentation)
import glob
from utils.utils import FileNameRetriever
import subprocess

                            
fnr = FileNameRetriever(bids_root)

def run_command(command):
    subprocess.run(command, shell=True, capture_output=True, check=True)

if do_hippocampus_segmentation:
    for subj in subjects:  
        subsubj = "sub-" + subj
        subjects_dir = fnr.get_filename(subsubj, "subjects_dir")
        print(f"Now running hippocampal segmentation for subject: {subj}\nThis might take some time")
        command = "segmentHA_T1.sh sub-" + subj + " " + subjects_dir
        run_command(command)



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
use freesurfers asegstats2table for subjects (bash) and extract from there
visualize single subjects segmentation and statistics:
    right to left difference?
    hippocampi too small relative to brain volume?
visualize single subject and statistics vs. grand average of healthy subjects in same age bin

--> doesn't work yet.

"""