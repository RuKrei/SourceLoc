#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)



import os
import shutil
import argparse
from os.path import join as opj
from dicom2nifti import convert_directory
import glob
import mne
from mne_bids import make_dataset_description, \
                        BIDSPath, write_anat, write_raw_bids
from utils import utils as u
import platform
from nipype.interfaces.freesurfer import ReconAll


# configuration
bids_root = "C:\\Users\\rudik\\MEG\\playground\\BIDS_root"
extras_directory = "C:\\Users\\rudik\\MEG\\playground\\extras"
input_folder = "C:\\Users\\rudik\\MEG\\playground\\input_folder"
openmp = n_jobs = 4
splitter = "\\" if platform.system().lower().startswith("win") else "/"
FS_SUBJECTS_DIR = os.environ.get("SUBJECTS_DIR")
if FS_SUBJECTS_DIR == None:
    print(f"It seems freesurfer is not properly set up on your computer")
    FS_SUBJECTS_DIR = "\\\\wsl.localhost\\Ubuntu-20.04\\usr\\local\\freesurfer\\7-dev\\subjects"
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject", action="store", 
                        type=str, required=False,
                        help="Name of the Patient/ Subject to process")
    parser.add_argument("--bidsroot", action="store", type=str, required=False, 
                        help="Specify BIDS root directory to use")
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
    #parser.add_argument("--openmp", action="store", type=str, required=False, 
    #                    help="Specify how many jobs/ processor cores to use")
    #parser.add_argument("--openmp", action="store", type=str, required=False, 
    #                    help="Specify how many jobs/ processor cores to use")
    args = parser.parse_args()  

# define subject
    subject = args.subject    
    if not subject:
        poss = [s for s in os.listdir(input_folder)]
        print(f"No subject specified, maybe you want to choose from those:\n {poss}")
        subject = input()
    
    print(f"Subject = {subject}")
    
# additional arguments
    if args.openmp:
        openmp = args.openmp
        print(f"Using {openmp} processor cores/ jobs.")
    
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


# make sure we have an MRI, convert, if necessary
    anafolder = opj(input_folder, subject)
    nii = glob.glob(opj(input_folder, subject, "*.nii*"))
    if os.path.isdir(anafolder) and nii == []:
        folder = str(glob.glob((anafolder + '/1*/100*/100*'), recursive=True)[0])
        try:
            convert_directory(
                dicom_directory=folder, 
                output_folder=anafolder, 
                compression=True, 
                reorient=True)
            print("\nMRI was converted to .nii.gz\n")
        except Exception as e:
            print(f"Something went wrong trying to convert the MRI to nifti: {e}")
    nii = glob.glob(opj(input_folder, subject, "*.nii*"))
    if nii == []:
        print("No anatomical data found, did you provide an MRI?")
    else:
        print(f"\nMRI already is in nifti-Format: {nii}\nDoing nothing...\n")

# Check if only freesurfer segmentation was desired and comply, if true
# Naturally, this only works with a freesurfer environment 
    if args.fsonly and args.fsonly.lower() == "true":
        command = f"recon-all -i {nii[0]} -s {subject} -openmp {openmp} - all"
        u.run_shell_command(command)
    
# create folder structure and copy 
    fbase = opj(bids_root, "derivatives", "sub-" + subject)
    fsrc = opj(fbase, "source_model")
    fanat = opj(fbase, "freesurfer")
    fprep = opj(fbase, "preprocessing")
    spikes = opj(fbase, "spikes")
    freq = opj(fbase, "frequency_distribution")
    conn = opj(fbase, "connectivity")
    ftrans = opj(fbase, "trans_files")
    freport = opj(fbase, "report")
    disclaimer = opj(extras_directory, "MEG_disclaimer.png")
    title = opj(extras_directory, "MEG_title.png")
    report_files = [disclaimer, title]
    ana_dir = opj(extras_directory, "anatomy_templates")
    ana_files = glob.glob(ana_dir)  
    # create folders
    folder_list = [fbase, fsrc, fanat, fprep, spikes, freq, conn, ftrans, freport]
    for fld in folder_list:
        if not os.path.exists(fld):
            os.makedirs(fld, exist_ok=True)
            print(f"Folder {fld} created.")
    # copy disclaimer and title for report
    for f in report_files:
        try:
            fi = f.split(splitter)[-1]
            target = os.path.join(freport, fi)
            shutil.copyfile(f, target)
        except Exception as e:
            print(e)
        #copy fsaverage + fsaverage_sym to local subjects anatomy folder
        for f in ana_files:
            if not (f.split(splitter)[-1].endswith(".png")):
                try:
                    u.recursive_overwrite(f, fanat)
                    print(f"Copying: {f}")
                except Exception as e:
                    print(e)
    
# create BIDS dataset
    raws = glob.glob(input_folder + "/*.fif")
    raws = [f for f in raws if subject in f]
    print(f"The following raw files were found:\n{raws}")
    bids_path = BIDSPath(subject=subject, session="resting", task="resting", 
                           root=bids_root, processing=None)
    fbase = os.path.join(bids_root, "derivatives", "sub-" + subject)
    fmeg = opj(fbase, "meg")

    # anatomy
    try:
        for n in nii:
            bids_path.update(root=bids_root)
            write_anat(n, bids_path=bids_path, overwrite=True)
    except Exception as e:
        print(e)

    # MEG
    # The following processing flags are added:
    # - .fif                                --> None
    # - tsss-trans.fif without event-file   --> tsssTrans
    # - tsss-trans.fif + event-file         --> tsssTransEve
    prepper = u.RawPreprocessor()
    # files are processed (Spike-Selection, Maxgfilter), so they should not 
    # be stored at BIDS-root anymore
    derivatives_root = opj(bids_root, "derivatives")  
    for run, rawfile in enumerate(raws):
        run +=1
        if "tsss" in rawfile:
            # --> search for matching eventfile and combine
            eve_name = rawfile.split(".fif")[0] + "_Events.csv"
            if not os.path.isfile(eve_name):
                eve_name = rawfile.split(".fif")[0] + "_Events.txt"
            if os.path.isfile(eve_name): # if fif-file matches event-file --> add events to fif-file
                print(f"\n\nNow adding Events ({eve_name}) to fif ({rawfile})\n\n")
                raw = mne.io.read_raw(rawfile, preload=True)
                event_file, event_dict = prepper.transform_eventfile(eve_name)
                raw.add_events(event_file)
                bids_path.update(root=derivatives_root, processing="tsssTransEve", run=run)
                raw.save(rawfile, overwrite=True)
                raw = mne.io.read_raw(rawfile, preload=False)
                write_raw_bids(raw, bids_path, event_id=event_dict, events_data=event_file.to_numpy(), overwrite=True)
            else: # tsss-file, but no events
                print(f"\n\nFound tsss-file {rawfile}, but no matching eventfile.\n\n")
                raw = mne.io.read_raw(rawfile, preload=False)
                bids_path.update(root=derivatives_root, processing="tsssTrans")
                write_raw_bids(raw, bids_path, overwrite=True)
        else: # unprocessed raw file
            print("\n\nFound raw file {rawfile}, saving in BIDS base.\n\n")
            bids_path.update(root=bids_root, processing=None)
            # write raw BIDS, as below
            raw = mne.io.read_raw(rawfile)
            write_raw_bids(raw, bids_path, overwrite=True)
    
    # Dataset
    the_roots = [bids_root, derivatives_root]
    for r in the_roots:
        make_dataset_description(r, 
                            name="CDK Epilepsy Dataset", 
                            data_license="closed", 
                            authors="Rudi Kreidenhuber", 
                            overwrite=True)

# Adjust subject variable, so we comply with BIDS convention
    if not subject.startswith("sub-"):
        subject = "sub-" + subject

# Freesurfer segmentation - nipype makes this very convenient
    fs_folder = opj(FS_SUBJECTS_DIR, subject)
    if os.path.isdir(fs_folder):
        print(f"A freesurfer segmentation for {subject} exists at:\n \
                {fs_folder}")
        
        
    else:
        anafolder = opj(bids_root, subject, "ses-resting", "anat")
        print(f"anafolder {anafolder}")
        nii = glob.glob(anafolder + splitter +  "*.nii*")[0]
        reconall = ReconAll()
        reconall.inputs.subject_id = subject
        reconall.inputs.T1_files = nii
        reconall.inputs.directive = 'all'
        reconall.inputs.subjects_dir = FS_SUBJECTS_DIR
        reconall.inputs.openmp = openmp
        reconall.inputs.flags = "-3T"
        reconall.run()
    this_subjects_dir = opj(fanat, subject)
    if not os.path.isdir(this_subjects_dir):
        u.recursive_overwrite(fs_folder, this_subjects_dir)
    else:
        print(f"Freesurfer segmentation already exists at {this_subjects_dir} \
                --> doing nothing...")
    
# Head model --> fails without freesurfer
    if not os.path.isfile(opj(this_subjects_dir, "bem", subject + "-head.fif")):
        try:    
            mne.bem.make_watershed_bem(subject=subject, subjects_dir=this_subjects_dir, overwrite=True)
        except Exception as e:
            print("#"*30)
            print("#"*30)
            print("#"*30)
            print(f"Failed to make watershed BEM for {subject}\n--> {e}")

# Cortex source space
    srcfilename = opj(fsrc, subject + "-" + spacing + "-src.fif")
    if not os.path.isfile(srcfilename):
        try:
            src = mne.setup_source_space(subject, spacing = spacing, 
                                            subjects_dir = fanat, 
                                            n_jobs=n_jobs, 
                                            verbose=True)
            mne.write_source_spaces(srcfilename, src, overwrite=True, verbose=True)
        except Exception as e:
            print("#"*30)
            print("#"*30)
            print("#"*30)
            print(f"Failed to setup source space with spacing {spacing} \
                    for {subject} --> {e}")


# Volume source space
    srcfilename = opj(fsrc, subject  + "-vol-src.fif")
    if not os.path.isfile(srcfilename):
        try:    
            src_vol = mne.setup_volume_source_space(subject, pos=3.0, 
                                            subjects_dir = fanat, 
                                            verbose=True)
            mne.write_source_spaces(srcfilename, src_vol, overwrite=True, verbose=True)
        except Exception as e:
            print("#"*30)
            print("#"*30)
            print("#"*30)
            print(f"Failed to setup Volume source space for {subject} --> {e}")







if __name__ == '__main__':
    main()


"""
To do:
    update auto-recon-all to also take care of head-model (--> watershed uses freesurfer)
    
"""