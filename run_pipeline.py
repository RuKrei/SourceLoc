#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)



import os
import argparse
from os.path import join as opj
import glob
import mne
from mne_bids import make_dataset_description, \
                        BIDSPath, write_anat, write_raw_bids, \
                        read_raw_bids
from src.utils import utils as u
from src import Anatomist, Folderer
import platform
import pickle
import matplotlib.pyplot as plt
import numpy as np


####################################################################
####################################################################
####################################################################
# configuration
#bids_root = "C:\\Users\\rudik\\MEG\\playground\\BIDS_root"
#extras_directory = "C:\\Users\\rudik\\MEG\\playground\\extras"
#input_folder = "C:\\Users\\rudik\\MEG\\playground\\input_folder"

bids_root = "/home/idrael/playground/BIDS_root"
extras_directory = "/home/idrael/playground/extras"
input_folder = "/home/idrael/playground/input_folder"

openmp = n_jobs = 8

# Filter and resample
l_freq: float = 0.1    # lower pass-band edge
h_freq: float = 70.   # higher pass-band edge
fir_design: str = "firwin"
s_freq: int = 300
# ECG artifact correction
ecg_channel = "ECG003"
n_grad: int =1
n_mag: int = 1
n_eeg: int = 1

# Frequency spectrum
freq_bands = dict(                #the frequency bands of interest for the analysis
                delta=(1, 4), 
                theta=(4, 7), 
                alpha=(8, 12), 
                beta=(13, 29), 
                gamma=(30, 80))

# Source localization
source_loc_methods = ["dSPM", "eLORETA"]
snr = 3
lambda2 = 1. / snr ** 2
peaks_tmin = -.025
peaks_tmax = 0.
peaks_mode = "abs"                                      # How to deal with the sign of the data:  "pos", "neg", or "abs"
peaks_nr_of_points = 5                                   

dip_times = {   'min20ms':  (-0.025,-0.020),            # Time points in Miliseconds for Equivalent current dipole fit
                'min15ms':  (-0.019,-0.015),
                'min10ms':  (-0.014,-0.010),
                'min5ms':  (-0.009,-0.005),
                'peak':     (-0.004,0.000)}

####################################################################
####################################################################
####################################################################



def main():
    splitter = "\\" if platform.system().lower().startswith("win") else "/"
    FS_SUBJECTS_DIR = os.environ.get("SUBJECTS_DIR")
    if FS_SUBJECTS_DIR == None:
        print(f"It seems freesurfer is not properly set up on your computer")
        # might also mean we are on windows, so:
        FS_SUBJECTS_DIR = "\\\\wsl.localhost\\Ubuntu-20.04\\usr\\local\\freesurfer\\7-dev\\subjects"
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject", action="store", 
                        type=str, required=False,
                        help="Name of the Patient/ Subject to process")
    parser.add_argument("--bidsroot", action="store", type=str, required=False, 
                        help="Specify a different BIDS root directory to use")
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

# define subject
    subject = args.subject    
    if not subject:
        poss = [s for s in os.listdir(input_folder)]
        print(f"No subject specified, maybe you want to choose from those:\n {poss}")
        subject = input()

    if not subject.startswith("sub-"):
        ject = str(subject)
        subject = "sub-" + subject
    else:
        ject = subject.split("sub-")[-1]
    print(f"Subject = {ject}")
    
# additional arguments
    if args.openmp:
        n_jobs = openmp = int(args.openmp)
        print(f"Using {openmp} processor cores/ jobs.")
    else:
        n_jobs = openmp = int(1)
    
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


# MRI to nii.gz, then freesurfer, then hippocampal subfields
# Naturally, this only works with a freesurfer environment 
# and this will take some time...
    anafolder = opj(input_folder, ject)
    if os.path.isdir(anafolder):
        rap = Anatomist.RawAnatomyProcessor(anafolder, FS_SUBJECTS_DIR, n_jobs=n_jobs)
        try:
            rap.run_anatomy_pipeline()
        except Exception as e:
            print(f"Something went wrong while processing anatomy: {e}")

# Check if only freesurfer segmentation was desired and comply, if true
    if args.fsonly and args.fsonly.lower() == "true":
        exit()
    
# create folder structure and copy 
    dfc = Folderer.DerivativesFoldersCreator(BIDS_root=bids_root, 
                                            extras_directory=extras_directory, 
                                            subject=subject)
    dfc.make_derivatives_folders()

# copy freesurfer files to local subjects_dir
    try:
        segmentation = opj(FS_SUBJECTS_DIR, subject)
        target = opj(dfc.fanat, subject)
        if not os.path.isdir(target):
            os.mkdir(target)
        print(f"Copying freesurfer segmentation {segmentation} to {target}")
        dfc._recursive_overwrite(segmentation, target)
    except Exception as e:
        print(e)

# create source models
    sourcerer = Anatomist.SourceModeler(subjects_dir=dfc.fanat, subject=subject, spacing=spacing, n_jobs=n_jobs)
    sourcerer.calculate_source_models()















# To do:
# Load MEG files, filter, resample, supress artifacts, create epochs and concatenate epochs
# save as -epo.fif
# MEG files are expected to be maxfiltered + transpositioned
# --> _tsss_trans.fif

    # parse list of appropriate raws
    raws = glob.glob(input_folder + "/*.fif")
    raws = [f for f in raws if ject in f]
    print(f"The following raw files were found:\n{raws}")
    prepper = u.RawPreprocessor()
    for run, rawfile in enumerate(raws):
        if "tsss" in rawfile and ject in rawfile:
            # --> search for matching eventfile and combine
                rawname = rawfile.strip(".fif") + "_prep.fif"
                if not os.path.isfile(rawname) and not "_prep" in rawfile:  # --> if this has not been done already
                    try:
                        raw, event_file, event_dict = prepper.combine_raw_and_eve(rawfile, run=run)  
                        event_file = event_file.to_numpy()  
                    except:  # this fails if no Event-file exists, or if it is corrupted
                        print(f"\n\nFound tsss-file {rawfile}, but no matching eventfile.\n\n")
                        event_file = None
                        event_dict = None
                        raw = mne.io.read_raw(rawfile, preload=False, on_split_missing="ignore")
                    # preprocessing
                    raw = prepper.filter_raw(raw, l_freq=l_freq, fir_design=fir_design,
                                                h_freq=h_freq, n_jobs=n_jobs)
                    raw = prepper.resample_raw(raw, s_freq=s_freq, n_jobs=n_jobs)

                    
                    
                    # Artifacts        
                    # ECG artifacts
                    # It's smarter to supervise this step (--> look at the topomaps!)
                    raw.load_data()
                    try:
                        ecg_projs, _ = mne.preprocessing.compute_proj_ecg(raw, n_grad=n_grad, n_mag=n_mag, 
                                                                          n_eeg=n_eeg, reject=None)
                        # lets not do this now......
                        raw.add_proj(ecg_projs, remove_existing=False)
                        fig = mne.viz.plot_projs_topomap(ecg_projs, info=raw.info, show=False)
                        savename = os.path.join(dfc.fprep, "ECG_projs_Topomap.png")
                        fig.savefig(savename)
                    except Exception as e:
                        print(e)
                        print("ECG - Atrifact correction failed!")
                    #EOG artifacts    
                    # It's a bad idea to do this in an automated step
                    try:
                        eog_evoked = mne.preprocessing.create_eog_epochs(raw).average()
                        #eog_evoked.apply_baseline((None, None))
                        eog_projs, _ = mne.preprocessing.compute_proj_eog(raw, n_grad=n_grad, n_mag=n_mag, n_eeg=n_eeg, 
                                                                    n_jobs=n_jobs, reject=None)
                        raw.add_proj(eog_projs, remove_existing=False) # --> don't do this in the early stages - see documentation
                        figs = eog_evoked.plot_joint(show=False)
                        for idx, fig in enumerate(figs):
                            savename = os.path.join(dfc.fprep, "EOG Topomap_" + str(idx) + ".png")
                            fig.savefig(savename)
                    except Exception as e:
                        print(e)
                        print("EOG - Atrifact correction failed!")
                    raw.apply_proj()
                    
                    
                    
                    
                    
                    
                    # save raw
                    raw.save(rawname, overwrite=True)
                    del(raw)
    # concatenate filtered and resampled files
    raws = glob.glob(input_folder + "/*.fif")
    raws = [f for f in raws if ject in f]
    raws = [f for f in raws if "_prep" in f]
    all_raws = dict()
    concatname = opj(os.path.dirname(raws[0]), str(subject) + "_concat.fif")
    if not os.path.isfile(concatname):
        for r in raws:
            all_raws[r] = mne.io.read_raw(r, preload=False)
        print(f"\n\n\n\nConcatenating files: {raws}")
        try:
            #raw = mne.concatenate_raws(list(all_raws[k] for k in all_raws.keys()))
            raw = mne.concatenate_raws([all_raws[r] for r in all_raws.keys()])
            print("Rawfiles have been concatenated....")

        except Exception as e:
            print(f"Failed trying to concatenate raw file\n {r} --> {e}")
            #print("Loading only first raw file!")
            #raw = mne.io.read_raw(raws[0])
        concatname = opj(os.path.dirname(raws[0]), str(subject) + "_concat.fif")
        print(f"Saving concatenated rawfile as {concatname}")
        raw.save(concatname)

    # Save in BIDS format
    derivatives_root = opj(bids_root, "derivatives")
    # meg
    bids_path = BIDSPath(subject=ject, session="resting", task="resting", 
                           root=derivatives_root, processing="concat")
    raw = mne.io.read_raw(concatname, preload=False)
    write_raw_bids(raw, bids_path, overwrite=True)
    # anatomy  
#    fbase = os.path.join(bids_root, "derivatives", "sub-" + ject)
    nii = glob.glob(opj(input_folder, ject, "*.nii*"))
    bids_path.update(root=bids_root)
    try:
        for n in nii:      
            write_anat(n, bids_path=bids_path, overwrite=True)
    except Exception as e:
        print(e)
    bids_path.update(root=bids_root)


    # Create Dataset
    the_roots = [bids_root, derivatives_root]
    for r in the_roots:
        make_dataset_description(r, 
                            name="CDK Epilepsy Dataset", 
                            data_license="closed", 
                            authors="Rudi Kreidenhuber", 
                            overwrite=True)


    """
    
    # files are processed (Spike-Selection, Maxgfilter), so they should not 
    # be stored at BIDS-root anymore
    all_raws = dict()
    for run, rawfile in enumerate(raws):
        run +=1
        if "tsss" in rawfile and ject in rawfile:
            # --> search for matching eventfile and combine
                try:
                    raw, event_file, event_dict = prepper.combine_raw_and_eve(rawfile, run=run)  
                    event_file = event_file.to_numpy()  
                except:  # this fails if no Event-file exists, or if it is corrupted
                    print(f"\n\nFound tsss-file {rawfile}, but no matching eventfile.\n\n")
                    event_file = None
                    event_dict = None
                    raw = mne.io.read_raw(rawfile, preload=True)
                #raw = prepper.filter_raw(raw, l_freq=l_freq, fir_design=fir_design,
                #                            h_freq=h_freq, n_jobs=n_jobs)
                raw = prepper.resample_raw(raw, s_freq=s_freq, n_jobs=n_jobs)
                all_raws[run] = raw.copy()
                bids_path.update(root=derivatives_root, processing="Prepared", run=run)
                modfile = rawfile.split(".fif")[0] + "-prep.fif"
                raw.save(modfile, overwrite=True)
                raw = mne.io.read_raw(modfile, preload=False)
                write_raw_bids(raw, bids_path, event_id=event_dict, events_data=event_file, overwrite=True)


    bids_derivatives = BIDSPath(subject=ject, 
                    datatype="meg", 
                    session="resting", 
                    task="resting", 
                    root=derivatives_root, 
                    processing="Concat", 
                    suffix="meg")
    raw.save(modfile, overwrite=True)
    raw = mne.io.read_raw(modfile, preload=False)
    write_raw_bids(raw, bids_derivatives, overwrite=True)
    del all_raws


#######################################################################################
#######################################################################################
#######################################################################################
# -->   Here would be the place to filter + concatenate input files and save one file
#       to be processed further
#######################################################################################
#######################################################################################
#######################################################################################
    





# .fif data preparations
# --> To do: concatenate, if multiple files
# Be less confusing with file loading logic
    # load all files
    raw = False
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
    if not raw:
        meg_dir = fbase + splitter + "ses-resting" + splitter + "meg"
        concat_fif = meg_dir + splitter + "*Concat*.fif"
        all_raws = glob.glob(concat_fif)
        if all_raws == []:
            tsssTrans_fif = meg_dir + splitter + "*tsssTrans*.fif"
            all_raws = glob.glob(tsssTrans_fif)
        print(f"\n\nThe following files are being processed: {all_raws}\n\n")
    
    # Check, if this has already been done
    meg_dir = fbase + splitter + "ses-resting" + splitter + "meg"
    eve_fif = meg_dir + splitter + "*tsssTransEve*.fif"
    eve_fif = glob.glob(eve_fif)
    noeve_fif = meg_dir + splitter + "*tsssTransNo*.fif"
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
#
        # filter
        #raw = prepper.filter_raw(raw, l_freq=l_freq, h_freq=h_freq, n_jobs=n_jobs)
        # resample
        print(f"Resample to {s_freq}")
        print(f"Original sampling frequency was: {raw.info['sfreq']}")
        raw = prepper.resample_raw(raw, s_freq=s_freq, n_jobs=n_jobs)
        
        




    """



        
    """    
        
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
    """

# Coregistration --> this doesn't work with WSLg - from here on
# run on windows, if you are on a windows machine
    transfile = opj(dfc.ftrans, subject + "-trans.fif")
    if os.path.isfile(transfile):
        print(f"Skipping coregistration, because a transfile ({transfile}) already exists")
    else:
        print(f"\n\n\n--> Transfile should be called: {transfile}\n\n\n")
        try:
            mne.gui.coregistration(subject=subject, subjects_dir=dfc.fanat, inst=bids_path, advanced_rendering=False) # BIDS: inst=raw.filenames[0])
        except:
            print("failed with bids_derivatives folder")
            rawfile = opj(dfc.fbase, "ses-resting", "meg") + splitter + "*concat_meg.fif"
            print(f"Rawfile = {rawfile}")
            rawfile = glob.glob(rawfile)[0]
            mne.gui.coregistration(subject=subject, subjects_dir=dfc.fanat, inst=rawfile, advanced_rendering=False)





# frequency spectrum
    bem_sol = opj(dfc.fsrc, subject + "-3-layer-BEM-sol.fif")
    fwd_name = opj(dfc.fsrc, subject + "-fwd.fif")
    srcfilename = opj(dfc.fsrc, subject + "-" + spacing + "-src.fif")
    filebase = str(subject) + "_Freqs"
    all_stcs_filename = (filebase + '-stc-psd-MNE.pkl')
    all_stcs_filename = opj(dfc.freq, all_stcs_filename)
    sensor_psd_filename = (filebase + '-sensor-psd-MNE.pkl')
    sensor_psd_filename = opj(dfc.freq, sensor_psd_filename)
    if not os.path.isfile(all_stcs_filename) or not os.path.isfile(sensor_psd_filename):  # so this should run only on the first file..
        raw.load_data()
        if os.path.isfile(fwd_name):
            fwd = mne.read_forward_solution(fwd_name)
        else:    
            fwd = mne.make_forward_solution(raw.info, src=srcfilename, bem=bem_sol,
                                            trans=transfile, 
                                            meg=True, eeg=False, mindist=0.2, 
                                            ignore_ref=False, 
                                            n_jobs=n_jobs, verbose=True)
            mne.write_forward_solution(fwd_name, fwd)
        noise_cov = mne.compute_raw_covariance(raw, method="empirical", n_jobs=n_jobs)
        inv = mne.minimum_norm.make_inverse_operator(raw.info, forward=fwd, noise_cov=noise_cov, 
                                            loose="auto", depth=0.8)
        snr = 3.
        lambda2 = 1. / snr ** 2
        stc_psd, sensor_psd = mne.minimum_norm.compute_source_psd(raw, inv, lambda2=lambda2, 
                                                method='MNE', 
                                                fmin=1, fmax=45, n_fft=2048, n_jobs=n_jobs, 
                                                return_sensor=True, verbose=True)
        pickle.dump(stc_psd, open(all_stcs_filename, "wb"))
        pickle.dump(sensor_psd, open(sensor_psd_filename, "wb"))
    else:
        stc_psd = pickle.load(open(all_stcs_filename, "rb"))
        sensor_psd = pickle.load(open(sensor_psd_filename, "rb"))
    # Visualization
    topos = dict()
    stcs = dict()
    topo_norm = sensor_psd.data.sum(axis=1, keepdims=True)
    stc_norm = stc_psd.sum()
    for band, limits in freq_bands.items():         # normalize...
        data = sensor_psd.copy().crop(*limits).data.sum(axis=1, keepdims=True)
        topos[band] = mne.EvokedArray(100 * data / topo_norm, sensor_psd.info)
        stcs[band] = 100 * stc_psd.copy().crop(*limits).sum() / stc_norm.data
    brain = dict()
    x_hemi_freq = dict()
    mne.viz.set_3d_backend('pyvista')
    for band in freq_bands.keys():
        brain[band] = u.plot_freq_band_dors(stcs[band], band=band, subject=subject, 
                                            subjects_dir=dfc.fanat, filebase=filebase)
        freqfilename3d = (filebase + '_' + band + '_freq_topomap_3d_dors.png')
        freqfilename3d = os.path.join(dfc.freq, freqfilename3d)
        image = brain[band].save_image(freqfilename3d)
        brain_lh, brain_rh = u.plot_freq_band_lat(stcs[band], band=band, subject=subject, 
                                                    subjects_dir=dfc.fanat,
                                                    filebase=filebase)                                           
        freqfilename3d = (filebase + '_' + band + '_freq_topomap_3d_lat_lh.png')
        freqfilename3d = os.path.join(dfc.freq, freqfilename3d)
        image = brain_lh.save_image(freqfilename3d)
        freqfilename3d = (filebase + '_' + band + '_freq_topomap_3d_lat_rh.png')
        freqfilename3d = os.path.join(dfc.freq, freqfilename3d)
        image = brain_rh.save_image(freqfilename3d)
        brain_lh, brain_rh = u.plot_freq_band_med(stcs[band], band=band, subject=subject, 
                                                    subjects_dir=dfc.fanat,
                                                    filebase=filebase)                                           
        freqfilename3d = (filebase + '_' + band + '_freq_topomap_3d_med_lh.png')
        freqfilename3d = os.path.join(dfc.freq, freqfilename3d)
        image = brain_lh.save_image(freqfilename3d)
        freqfilename3d = (filebase + '_' + band + '_freq_topomap_3d_med_rh.png')
        freqfilename3d = os.path.join(dfc.freq, freqfilename3d)
        image = brain_rh.save_image(freqfilename3d)
    # 2. Cross hemisphere comparison
        # make sure fsaverage_sym exists in local subjects dir:
        target = os.path.join(dfc.fanat, "fsaverage_sym")
        if not os.path.isdir(target):
            # try to find it in $SUBJECTS_DIR and copy
            os_subj_dir = os.environ.get("SUBJECTS_DIR")
            fs_avg_sym_dir = os.path.join(os_subj_dir, "fsaverage_sym")
            u.recursive_overwrite(fs_avg_sym_dir, target)
        mstc = stcs[band].copy()
        mstc = mne.compute_source_morph(mstc, subject, 'fsaverage_sym',
                                                    smooth=5,
                                                    warn=False,
                                                    subjects_dir=dfc.fanat).apply(mstc)
        morph = mne.compute_source_morph(mstc, 'fsaverage_sym', 'fsaverage_sym',
                                                    spacing=mstc.vertices, warn=False,
                                                    subjects_dir=dfc.fanat, xhemi=True,
                                                    verbose='error')
        stc_xhemi = morph.apply(mstc)
        diff = mstc - stc_xhemi
        title = ('blue = RH; ' + subject + ' -Freq-x_hemi- ' + band)
        x_hemi_freq[band] = diff.plot(hemi='lh', subjects_dir=dfc.fanat, 
                                size=(1200, 800), time_label=title,
                                add_data_kwargs=dict(time_label_size=10))
        freqfilename3d = (filebase + '_x_hemi_' + band + '.png')
        freqfilename3d = os.path.join(dfc.freq, freqfilename3d)
        image = x_hemi_freq[band].save_image(freqfilename3d)






# Source localization
    #target_dir = os.path.join(derivatives_root, subject, 
    #                        "ses-resting", "meg", subject)
    #all_raws = glob.glob(target_dir + "*tsssTransEve_meg.fif")    # should already be concatenated
    #raw = read_raw_bids(all_raws[0])
    

    bids_derivatives.update(processing="finalEpochs", session="99")
    fname = os.path.join(derivatives_root, subject, "ses-resting", "meg", bids_derivatives.basename)
    fname = fname + ".fif"
    raw.save(fname, overwrite=True)
    events, event_ids = mne.events_from_annotations(raw)

    epochs = mne.Epochs(raw, events=events, event_id=event_ids, tmin=-1.5, tmax=1, 
                        baseline=(-1.5,-1), on_missing = "ignore",
                        event_repeated="merge")
    noise_cov_file = opj(spikes, "Spikes_noise_covariance.pkl")
    if not os.path.isfile(noise_cov_file):
        noise_cov = mne.compute_covariance(epochs, tmax=-1., 
                                    method='auto',
                                    n_jobs=n_jobs,
                                    rank="full")
        pickle.dump(noise_cov, open(noise_cov_file, "wb"))
    else:
        with open(noise_cov_file, 'rb') as f:
            noise_cov = pickle.load(f)
    
    #data_cov = mne.compute_covariance(epochs,
    #                                tmin=-0.5, 
    #                                tmax=0.3, 
    #                                method='auto',
    #                                n_jobs=n_jobs)
    
    print(f"event_ids.keys = {event_ids.keys()}")
    for event in event_ids.keys():
        eventname = str(event)
        if eventname == "ignore_me" or eventname == "AAA" or eventname == ".ungrouped":
            print(f"Omitting event {event}")
        else:
            try:
                print(f"\n\n\nNow processing event: {event}")
                e = epochs[eventname].load_data().average()
                e_folder = os.path.join(spikes, eventname)
                cp_folder = os.path.join(spikes, eventname, "custom_pics")
                cts_folder = os.path.join(spikes, eventname, "custom_time_series")
                gp_folder = os.path.join(spikes, eventname, "generic_pics")
                folders = [e_folder, cp_folder, cts_folder, gp_folder]
                if not os.path.isdir(e_folder):
                    for f in folders:
                        os.mkdir(f)
                src = mne.read_source_spaces(srcfilename)
                bem_sol = opj(fsrc, subject + "-3-layer-BEM-sol.fif")

                fwd = mne.make_forward_solution(e.info, src=src, bem=bem_sol,
                                            trans=transfile, 
                                            meg=True, eeg=False, mindist=0.2, 
                                            ignore_ref=False, 
                                            n_jobs=n_jobs, verbose=True)
                inv = mne.minimum_norm.make_inverse_operator(e.info, forward=fwd, noise_cov=noise_cov, 
                                            loose=0.2, depth=0.8)
                # Distributed source models
                for m in source_loc_methods:
                    stc_name = "stc_" + m + "_" + eventname
                    stc_name = 'stc_' + m
                    if m == 'dSPM':
                        stc_name = mne.minimum_norm.apply_inverse(e, inv, lambda2,
                                    method='dSPM', pick_ori='vector')
                        surfer_kwargs = dict(hemi='split', subjects_dir=fanat,
                                    clim=dict(kind='percent', lims=[90, 96, 99.85]),
                                    views=['lat', 'med'], 
                                    colorbar=True,
                                    initial_time=0, time_unit='ms', 
                                    size=(1000, 800), smoothing_steps=10)
                        brain = stc_name.plot(**surfer_kwargs)
                        label = str(subject + " - " + eventname + " - Vector solution")
                        brain.add_text(0.1, 0.9, label, 'title', font_size=10)
                        img_f_name = ('img_stc_' + subject + '_' + eventname + '_' + m + '.png')
                        img_f_name = os.path.join(gp_folder, img_f_name)
                        brain.save_image(img_f_name)
                        stc_f_name = ('stc_' + subject + '_' + eventname + '_' + m)
                        stc_f_name = os.path.join(e_folder, stc_f_name)
                        stc_name.save(stc_f_name)
                    else:
                        stc_name = mne.minimum_norm.apply_inverse(e, inv, lambda2,
                                    method=m, pick_ori=None)
                        surfer_kwargs = dict(hemi='split', subjects_dir=fanat,
                                    clim=dict(kind='percent', lims=[90, 96, 99.85]),
                                    views=['lat', 'med'], 
                                    colorbar=True,
                                    initial_time=0, time_unit='ms', 
                                    size=(1000, 800), smoothing_steps=10)
                        brain = stc_name.plot(**surfer_kwargs)
                        label = str(subject + " - " + eventname + " - " +  m)
                        brain.add_text(0.1, 0.9, label, 'title', font_size=10)
                        img_f_name = ('img_stc_' + subject + '_' + eventname + '_' + m + '.png')
                        img_f_name = os.path.join(gp_folder, img_f_name)
                        brain.save_image(img_f_name)
                        stc_f_name = ('stc_' + subject + '_' + eventname + '_' + m)
                        stc_f_name = os.path.join(e_folder, stc_f_name)
                        e.save(stc_f_name)
                        if m == "eLORETA":
                            rh_peaks = u.get_peak_points(stc_name, hemi='rh', tmin=peaks_tmin, 
                                                        tmax=peaks_tmax, nr_points=peaks_nr_of_points, mode=peaks_mode)
                            lh_peaks = u.get_peak_points(stc_name, hemi='lh', tmin=peaks_tmin, 
                                                        tmax=peaks_tmax, nr_points=peaks_nr_of_points, mode=peaks_mode)
                            label = str(subject + " - " + eventname + " - " +  m + " - max. activation points")
                            brain.add_text(0.1, 0.9, label, 'title', font_size=10)
                            for p in rh_peaks:
                                brain.add_foci(p, color='green', coords_as_verts=True, hemi='rh', scale_factor=0.6, alpha=0.9)
                            for p in lh_peaks:
                                brain.add_foci(p, color='green', coords_as_verts=True, hemi='lh', scale_factor=0.6, alpha=0.9)
                            stc_f_name = ('stc_' + subject + '_' + eventname + '_' + m + "_with_peaks")
                            stc_f_name = os.path.join(e_folder, stc_f_name)
                            stc_name.save(stc_f_name)
                            img_f_name = ('img_stc_' + subject + '_' + eventname + '_' + m + '_with_peaks.png')
                            img_f_name = os.path.join(gp_folder, img_f_name)
                            brain.save_image(img_f_name)
                # Dipoles
                for start, stop in dip_times.values():
                    dip_epoch = e.copy().crop(start, stop).pick('meg')
                    ecd = mne.fit_dipole(dip_epoch, noise_cov, bem_sol, trans=transfile)[0]
                    best_idx = np.argmax(ecd.gof)
                    best_time = ecd.times[best_idx]
                    trans = mne.read_trans(transfile)
                    mri_pos = mne.head_to_mri(ecd.pos, mri_head_t=trans, subject=subject, subjects_dir=fanat)
                    t1_file_name = os.path.join(fanat, subject, 'mri', 'T1.mgz')
                    stoptime = str(abs(int(stop*1000)))
                    if stoptime == "5":
                        stoptime = "05"
                    title = str(eventname + ' - ECD @ minus ' + stoptime + ' ms')
                    t1_fig = u.plot_anat(t1_file_name, cut_coords=mri_pos[0], title=title)
                    t1_f_name_pic = ('img_ecd_' + eventname + '_' + '_Dipol_' + stoptime + '.png')
                    t1_f_name_pic = os.path.join(e_folder, "generic_pics", t1_f_name_pic)
                    t1_fig.savefig(t1_f_name_pic)
                    fig_3d = ecd.plot_locations(trans, subject, fanat, mode="orthoview")
                    fig_3d_pic = ('img_3d_ecd_' + eventname + '_' + '_Dipol_' + stoptime + '.png')
                    fig_3d_pic = os.path.join(e_folder, "generic_pics", fig_3d_pic)
                    fig_3d.savefig(fig_3d_pic)
                    plt.close("all")
            except Exception as e:
                print(e)
  

if __name__ == '__main__':
    main()


"""
To do:
    update auto-recon-all to also take care of head-model (--> watershed uses freesurfer)
    
    resample and concatenate at the beginning of the pipeline --> check, that it works properly

    fix filtering, resampling
    
"""