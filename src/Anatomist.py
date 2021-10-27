#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from nipype.interfaces.freesurfer import ReconAll
from dicom2nifti import convert_directory
import os
from os.path import join as opj
import argparse
import glob
import subprocess
import shutil

class RawAnatomyProcessor:
    """ This class is responsible for processing 
        MRIs. 
        It is the host for helper-functions that transform: 
        - DICOM  to nifti
        - nifti to freesurfer (recon-all)
        It also runs the watershed algorithm needed for the head model.
    """
    def __init__(self, mri_folder, subjects_dir, n_jobs):
        self.mri_folder = mri_folder
        self.subjects_dir = subjects_dir
        self.subject = os.path.dirname(self.subjects_dir)
        self.FS_SUBJECTS_DIR = os.environ.get("SUBJECTS_DIR")
        self.n_jobs = n_jobs
    
    def _dicom_to_nii(self):
        folder = str(glob.glob(opj(self.mri_folder, "1*", "100*", "100*"), recursive=True)[0])
        try:
            convert_directory(
                dicom_directory=folder, 
                output_folder=self.mri_folder, 
                compression=True, 
                reorient=True)
            print("\nMRI was converted to .nii.gz\n")
        except Exception as e:
            print(f"Something went wrong trying to convert the MRI to nifti: {e}")

    def _nii_to_freesurfer(self):
        reconall = ReconAll()
        if not self.subject.startswith("sub-"):
            self.subject = "sub-" + self.subject
        # check if freesurfer segmentation was already performed
        freesurfered = os.path.join(self.FS_SUBJECTS_DIR, self.subject)
        if not os.path.isdir(freesurfered):
            nii_file = glob.glob(opj(self.mri_folder, "*.nii*"))[0]
            reconall.inputs.subject_id = self.subject
            reconall.inputs.T1_files = nii_file
            reconall.inputs.directive = 'all'
            reconall.inputs.subjects_dir = self.FS_SUBJECTS_DIR
            reconall.inputs.openmp = self.n_jobs
            reconall.inputs.flags = "-3T"
            reconall.run()
        else:
            print(f"A freesurfer segmentation of subject {self.subject} already exists in {self.FS_SUBJECTS_DIR}")

    def _run_shell_command(command):
        subprocess.run(command, shell=True, capture_output=True, check=True)
    
    def _segment_hippocampal_subfields(self):
        """Only works with Matlab (or Matlab runtime environment)
        """
        hippofile = opj(self.FS_SUBJECTS_DIR, self.subject, "mri", "lh.hippoSfVolumes-T1.v21.txt")
        if not os.path.isfile(hippofile):
            try:
                print(f"Now running hippocampal segmentation for subject: {self.subject}\nThis might take some time")
                hipposeg = "segmentHA_T1.sh " + self.subject
                self._run_shell_command(hipposeg)
            except Exception as e:
                print(f"Something went wrong with the hippocampal segmentation! --> {e}")
        else:
            print(f"Omitting hippocampal segmentation for subject {self.subject}, as it already exists")

    def run_anatomy_pipeline(self):
        self._dicom_to_nii()
        self._nii_to_freesurfer()
        self._segment_hippocampal_subfields()
      

class HeadModeler:
    """ Generates desired head models.
        Expects to find subject of interest as processed
        with freesurfer in $SUBJECTS_DIR
    """
    pass

class HippocampusSegmenter:
    """Adds hippocampal subfield segmentation
    """
    pass


def main():
    FS_SUBJECTS_DIR = os.environ.get("SUBJECTS_DIR")
    if FS_SUBJECTS_DIR == None:
        print(f"It seems freesurfer is not properly set up on your computer")
    parser = argparse.ArgumentParser()    
    parser.add_argument("--subject", action="store", 
                        type=str, required=False,
                        help="Name of the Patient/ Subject to process")    
    parser.add_argument("--inputfolder", action="store", type=str, required=False, 
                        help="Specify a different data input folder")
    parser.add_argument("--fsonly", action="store", type=str, required=False, 
                        help="Use --fsonly true if you only want to do a freesurfer segmentation")      # do only freesurfer segmentation
    parser.add_argument("--openmp", action="store", type=str, required=False, 
                        help="Specify how many jobs/ processor cores to use")
    parser.add_argument("--srcspacing", action="store", type=str, required=False, 
                        help="Source spacing: \
                            -defaults to ico4 --> 2562 Source points \
                            || other options: \
                            oct5 --> 1026 Source points \
                            || oct6 --> 4098 Source points \
                            || ico5 --> 10242 Source points")    
    
    args = parser.parse_args() 
    
    subject = args.subject    
    if not subject:
        print(f"No valid subject specified, aborting!")
    print(f"Subject = {subject}")
    
    if args.openmp:
        openmp = n_jobs = args.openmp
        print(f"Using {openmp} processor cores/ jobs.")
    else:
        openmp = n_jobs = args.openmp
    
    if not args.srcspacing:
        spacing = "ico4"
    else:
        spacing = args.srcspacing
        if spacing in ["ico4", "oct5", "oct6", "ico5"]:
            print(f"Desired source spacing is {spacing}")
        else:
            print('The desired spacing isn\'t allowed, typo?\n \
                Options are: "ico4", "oct5", "oct6", "ico5"')
            raise Exception
    
    rap = RawAnatomyProcessor(subject, FS_SUBJECTS_DIR, n_jobs=1)
    rap.run_anatomy_pipeline()




if __name__ == '__main__':
    main()
