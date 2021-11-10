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
import mne

class RawAnatomyProcessor:
    """ This class is responsible for processing 
        MRIs. 
        It is the host for helper-functions that transform: 
        - DICOM  to nifti
        - nifti to freesurfer (recon-all)
        It also runs the watershed algorithm needed for the head model.
    """
    def __init__(self, mri_folder, FS_SUBJECTS_DIR, n_jobs):
        self.mri_folder = os.path.expanduser(mri_folder)
        self.subject = os.path.basename(os.path.normpath(self.mri_folder))
        if not self.subject.startswith("sub-"):
            self.subject = "sub-" + self.subject
        self.FS_SUBJECTS_DIR = FS_SUBJECTS_DIR
        self.n_jobs = int(n_jobs)
    
    def _dicom_to_nii(self):
        folder = opj(self.mri_folder, "1*", "100*", "1*")
        print(f"dicom_to_nii - folder = {folder}")
        folder = str(glob.glob(folder, recursive=True)[0])
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

    def _run_shell_command(self, command):
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
    
    def _run_watershed(self):
        try:
            mne.bem.make_watershed_bem(self.subject, self.FS_SUBJECTS_DIR, copy=True)
        except Exception as e:
            print(e)

    def run_anatomy_pipeline(self):
        self._dicom_to_nii()
        self._nii_to_freesurfer()
        self._segment_hippocampal_subfields()
        self._run_watershed()
      

class SourceModeler:
    """ Generates desired head models.
        Expects to find subject of interest as processed
        with freesurfer in $SUBJECTS_DIR
    """
    def __init__(self, subjects_dir=None, subject=None, spacing="ico4", n_jobs=1):
        self.subjects_dir = subjects_dir
        self.subject = subject
        self.fsrc = opj(self.subjects_dir, "..", "source_model")
        self.spacing = spacing
        self.n_jobs = n_jobs
    
    def _make_vol_source_space(self):
        srcfilename = opj(self.fsrc, self.subject  + "-vol-src.fif")
        if not os.path.isfile(srcfilename):
            try:    
                src_vol = mne.setup_volume_source_space(self.subject, pos=3.0, 
                                                subjects_dir = self.subjects_dir, 
                                                verbose=True)
                mne.write_source_spaces(srcfilename, src_vol, overwrite=True, verbose=True)
            except Exception as e:
                print(f"Failed to setup volume source space for {self.subject} --> {e}")
    
    def _make_cortex_source_space(self):
        srcfilename = opj(self.fsrc, self.subject  + "-" + self.spacing + "-src.fif")
        if not os.path.isfile(srcfilename):
            try:
                src = mne.setup_source_space(self.subject, spacing = self.spacing, 
                                                subjects_dir = self.subjects_dir, 
                                                n_jobs=self.n_jobs, 
                                                verbose=True)
                mne.write_source_spaces(srcfilename, src, overwrite=True, verbose=True)
            except Exception as e:
                print(f"Failed to setup source space with spacing {spacing} \
                        for {self.subject} --> {e}")


    def _make_bem_solution(self):
        bem_save_name = opj(self.fsrc, self.subject + "-single-shell-model")
        if not os.path.isfile(bem_save_name):
            try:
                bem = mne.make_bem_model(self.subject, ico=4, 
                                conductivity=[0.3],   
                                subjects_dir=self.subjects_dir, verbose=True)
                mne.write_bem_surfaces(bem_save_name, bem, overwrite=True) #, overwrite=True)
            except Exception as e:
                print(f"Failed to setup single shell BEM model for {self.subject} --> {e}")
        bem_sol_filename = opj(self.fsrc, self.subject + "-single-shell-BEM-sol.fif")
        if not os.path.isfile(bem_sol_filename):
            try:
                bem = mne.read_bem_surfaces(bem_save_name)    
                bem_sol = mne.make_bem_solution(bem)
                mne.write_bem_solution(bem_sol_filename, bem_sol, overwrite=True)
            except Exception as e:
                print(f"Failed to calculate BEM solution (single-shell) for {self.subject} --> {e}")

        # BEM Solutions - 3-layer-BEM
        bem_save_name = opj(self.fsrc, self.subject + "-3-layer-BEM-model.fif")
        if not os.path.isfile(bem_save_name):
            try:
                bem = mne.make_bem_model(self.subject, ico=4, 
                                conductivity=[0.3, 0.006, 0.3],   
                                subjects_dir=self.subjects_dir, verbose=True)
                mne.write_bem_surfaces(bem_save_name, bem, overwrite=True) #, overwrite=True)
            except Exception as e:
                print(f"Failed to calculate 3-layer BEM model for {self.subject} --> {e}")
        bem_sol_filename = opj(self.fsrc, self.subject + "-3-layer-BEM-sol.fif")
        if not os.path.isfile(bem_sol_filename):
            try:
                bem = mne.read_bem_surfaces(bem_save_name)    
                bem_sol = mne.make_bem_solution(bem)
                mne.write_bem_solution(bem_sol_filename, bem_sol, overwrite=True)
            except Exception as e:
                print(f"Failed to calculate 3-layer BEM solution for {self.subject} --> {e}")
                print("This is bad, please look into the freesurfer segmentation...")
                print("Alternatively, you might be able to run the analysis with a single-shell-head-model (look into the configuration file")

    def calculate_source_models(self):
        self._make_vol_source_space()
        self._make_cortex_source_space()
        self._make_bem_solution()


def main(inputfolder=None):
    """Transforms: 
        - DICOM  to nifti
        - nifti to freesurfer (recon-all)
        - runs watershed algorithm needed for head modeling.

    Args:
        inputfolder ([Directory path], mandatory): [Path to dicom or nii directory]. Defaults to None.

    Returns:
        None
        .nii.gz - file is stored in <inputfolder>
        freesurfer output is stored in $SUBJECTS_DIR
    """
    parser = argparse.ArgumentParser()     
    parser.add_argument("--inputfolder", action="store", type=str, required=True, 
                        help="Specify the mri data input folder")
    parser.add_argument("--SUBJECTS_DIR", action="store", type=str, required=False, 
                        help="Freesurfer subjects dir ($SUBJECTS_DIR)")
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
    
    inputfolder = args.inputfolder
    print(f"Inputfolder = {inputfolder}")
    subject = os.path.basename(os.path.normpath(inputfolder))
    print(f"Subjectname = {subject}")
    
    if args.openmp:
        openmp = n_jobs = int(args.openmp)
        print(f"Using {openmp} processor cores/ jobs.")
    else:
        openmp = n_jobs = int(1)

    if args.SUBJECTS_DIR:
        FS_SUBJECTS_DIR = args.SUBJECTS_DIR
    else:
        FS_SUBJECTS_DIR = os.environ.get("SUBJECTS_DIR")
    
    if FS_SUBJECTS_DIR == None:
        print(f"It seems freesurfer is not properly set up on your computer")

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
    
    rap = RawAnatomyProcessor(inputfolder, FS_SUBJECTS_DIR, n_jobs=n_jobs)
    rap.run_anatomy_pipeline()




if __name__ == '__main__':
    main()
