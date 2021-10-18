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


# configuration
bids_root = "C:\\Users\\rudik\\MEG\\playground\\BIDS_root"
extras_directory = "C:\\Users\\rudik\\MEG\\playground\\extras"
input_folder = "C:\\Users\\rudik\\MEG\\playground\\input_folder"
openmp = n_jobs = 4
splitter = "\\" if platform.system().lower().startswith("win") else "/"


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

# make sure we have an MRI, convert, if necessary
    anafolder = opj(input_folder, subject)
    nii = glob.glob(opj(input_folder, subject, "*.nii.gz"))
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
    nii = glob.glob(opj(input_folder, subject, "*.nii.gz"))
    if nii == []:
        print("No anatomical data found, did you provide an MRI?")
    else:
        print(f"\nMRI already is in .nii.gz-Format: {nii}\nDoing nothing...\n")

# Check if only freesurfer segmentation was desired and comply, if true
# Naturally, this only works with a freesurfer environment 
    if args.fsonly.lower() == "true":
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








if __name__ == '__main__':
    main()


"""
To do:
    ignore MRI, if freesurfer subjetc already exists
    mri-only option
    
"""