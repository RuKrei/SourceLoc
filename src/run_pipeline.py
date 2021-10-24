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
from mne.filter import _filter_attenuation
from mne_bids import make_dataset_description, \
                        BIDSPath, write_anat, write_raw_bids, \
                        read_raw_bids
from utils import utils as u
import platform
from nipype.interfaces.freesurfer import ReconAll


# configuration
#bids_root = "C:\\Users\\rudik\\MEG\\playground\\BIDS_root"
#extras_directory = "C:\\Users\\rudik\\MEG\\playground\\extras"
#input_folder = "C:\\Users\\rudik\\MEG\\playground\\input_folder"

bids_root = "/home/idrael/playground/BIDS_root"
extras_directory = "/home/idrael/playground/extras"
input_folder = "/home/idrael/playground/input_folder"

openmp = n_jobs = 8
splitter = "\\" if platform.system().lower().startswith("win") else "/"
FS_SUBJECTS_DIR = os.environ.get("SUBJECTS_DIR")
# Filter and resample
l_freq: int = 1, 
h_freq: int = 70, 
fir_design = "firwin"
s_freq: int = 300
# ECG artifact correction
ecg_channel = "ECG003"
n_grad: int =1
n_mag: int = 1
n_eeg: int = 1


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

# .fif data preparations
    # load all files
    bids_derivatives = BIDSPath(subject=subject.split("sub-")[-1], 
                        datatype="meg", 
                        session="resting", 
                        task="resting", 
                        run="01",
                        root=derivatives_root,
                        suffix="meg",
                        extension="fif") 
    bids_derivatives = bids_derivatives.update(processing="tsssTransEve")
    all_raws = bids_derivatives.match() 
    print(f"all_raws according to bids: {all_raws}")
    try:
        raw = read_raw_bids(bids_derivatives)
    except Exception as e:
        print(f"Couldn't load BIDS fif with events: {e}")
        bids_derivatives = bids_derivatives.update(processing="tsssTrans")
        all_raws = bids_derivatives.match() 
        try:
            raw = read_raw_bids(bids_derivatives)
        except Exception as e:
            print(f"Couldn't load any BIDS fif: {e}")

    if all_raws == []:
        meg_dir = fbase + splitter + "ses-resting" + splitter + "meg"
        concat_fif = meg_dir + splitter + "*Concat*.fif"
        all_raws = glob.glob(concat_fif)
        if all_raws == []:
            tsssTrans_fif = meg_dir + splitter + "*tsssTrans*.fif"
            all_raws = glob.glob(tsssTrans_fif)
        print(f"\n\nThe following files are being processed: {all_raws}\n\n")
    
    # Check, if this has already been done
    meg_dir = fbase + splitter + "ses-resting" + splitter + "meg"
    eve_fif = meg_dir + splitter + "*tsssTransEve_*.fif"
    eve_fif = glob.glob(eve_fif)
    noeve_fif = meg_dir + splitter + "*tsssTrans_*.fif"
    noeve_fif = glob.glob(noeve_fif)   
    if (eve_fif == [] and noeve_fif == []):
        # concatenate raws, if there is more than one - because why not.
        if (len(all_raws) > 1):
            try:
                raws = dict()
                for num, rawfile in enumerate(all_raws):
                    raws[num] = mne.io.read_raw(rawfile)
                raw = mne.concatenate_raws(list(raws.values()))
                print("Rawfiles have been concatenated....")
                bids_derivatives = BIDSPath(subject=subject.split("sub-")[-1], 
                                datatype="meg", 
                                session="resting", 
                                task="resting", 
                                root=derivatives_root, 
                                processing="tsssTransEveConcat", 
                                suffix="meg")
                write_raw_bids(raw, bids_derivatives, overwrite=True)

            except Exception as e:
                print(f"Failed trying to concatenate multiple raw files --> {e}")
                print("Loading only first raw file!")
                raw = mne.io.read_raw(all_raws[0])
        else:     
        # raw = read_raw_bids(bids_derivatives)   # this fails again, using bare MNE to load data file
            if not raw:
                raw = mne.io.read_raw(all_raws[0])

        # filter
        #raw = prepper.filter_raw(raw, l_freq=l_freq, h_freq=h_freq, n_jobs=n_jobs)
        # resample
        print(f"Resample to {s_freq}")
        print(f"Original sampling frequency was: {raw.info['sfreq']}")
        raw = prepper.resample_raw(raw, s_freq=s_freq, n_jobs=n_jobs)
        # ECG artifacts
        # It's smarter to supervise this step (--> look at the topomaps!)
        try:
            ecg_projs, _ = mne.preprocessing.compute_proj_ecg(raw, n_grad=n_grad, n_mag=n_mag, 
                                                              n_eeg=n_eeg, reject=None)
            # lets not do this now......
            raw.add_proj(ecg_projs, remove_existing=False)
            #raw.apply_proj(ecg_projs, verbose=None) # don't do this in the early stages - see documentation
            fig = mne.viz.plot_projs_topomap(ecg_projs, info=raw.info, show=False)
            savename = os.path.join(fprep, "ECG_projs_Topomap.png")
            fig.savefig(savename)
        except Exception as e:
            print(e)
            print("ECG - Atrifact correction failed!")
        #EOG artifacts    
            ## It's a bad idea to do this in an automated step
            #try:
            #    eog_evoked = mne.preprocessing.create_eog_epochs(raw).average()
            #    #eog_evoked.apply_baseline((None, None))
            #    eog_projs, _ = mne.preprocessing.compute_proj_eog(raw, n_grad=n_grad, n_mag=n_mag, n_eeg=n_eeg, 
            #                                                n_jobs=n_jobs, reject=None)
            #    raw.add_proj(eog_projs, remove_existing=False) # .apply_proj() --> don't do this in the early stages - see documentation
            #    figs = eog_evoked.plot_joint(show=False)
            #    for idx, fig in enumerate(figs):
            #        savename = os.path.join(preproc_folder, "EOG Topomap_" + str(idx) + ".png")
            #        fig.savefig(savename)
            #except Exception as e:
            #    print(e)
            #    print("EOG - Atrifact correction failed!")
        # save - if events have been found
        events, event_ids = mne.events_from_annotations(raw)

        bids_derivatives = BIDSPath(subject=subject.split("sub-")[-1], 
                                datatype="meg", 
                                session="resting", 
                                task="resting", 
                                root=derivatives_root, 
                                processing="tsssTransEvePreproc", 
                                suffix="meg")  
        if len(events) > 0:
            raw_temp = os.path.join(fprep, "temp.fif")
            raw.save(raw_temp, overwrite=True)
            raw = mne.io.read_raw(raw_temp, preload=False)                    
            write_raw_bids(raw, bids_derivatives, overwrite=True)
        else:
        # save - if no events
            raw_temp = os.path.join(fprep, "temp.fif")
            raw.save(raw_temp, overwrite=True)
            bids_derivatives.update(processing="tsssTransNoEvePreproc", run=run)   
            raw = mne.io.read_raw(raw_temp, preload=False)                    
            write_raw_bids(raw, bids_derivatives, overwrite=True)
    else:
        print("Omitting preprocessing steps, as preprocessed file has been found.")

# Coregistration --> this doesn't work with WSLg
    transfile = opj(ftrans, subject + "-trans.fif")
    if os.path.isfile(transfile):
        print(f"Skipping coregistration, because a transfile ({transfile}) already exists")
    else:
        print(f"\n\n\n--> Transfile should be called: {transfile}\n\n\n")
        try:
            mne.gui.coregistration(subject=subject, subjects_dir=fanat, inst=bids_derivatives, advanced_rendering=False) # BIDS: inst=raw.filenames[0])
        except:
            print("failed with bids_derivatives")
            rawfile = meg_dir + splitter + "*Preproc*.fif"
            print(f"Rawfile = {rawfile}")
            rawfile = glob.glob(rawfile)[0]
            mne.gui.coregistration(subject=subject, subjects_dir=fanat, inst=rawfile, advanced_rendering=False)




if __name__ == '__main__':
    main()


"""
To do:
    update auto-recon-all to also take care of head-model (--> watershed uses freesurfer)
    
    resample and concatenate at the beginning of the pipeline
    
"""