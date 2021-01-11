#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from configuration import (subjects, n_jobs, bids_root, data_root, session, 
                            concat_raws, l_freq, h_freq, fir_design, s_freq,
                            do_filter, do_resample, do_ecg_correction_ssp, 
                            do_ecg_correction_ica, do_ecg_correction_regression,
                            do_eog_correction_ssp, do_eog_correction_ica, 
                            do_eog_correction_regression, do_source_loc, 
                            n_grad, n_mag, n_eeg, eog_channel, ecg_channel,
                            do_source_loc)
import glob
import mne
from utils.utils import FileNameRetriever , RawPreprocessor
import utils.utils as u
import pandas as pd
from shutil import copyfile
from mne_bids import BIDSPath, write_raw_bids, read_raw_bids
import matplotlib.pyplot as plt

fnr = FileNameRetriever(bids_root)
prepper = RawPreprocessor()


for subj in subjects:
    if not subj.startswith("sub-"):
        subj = "sub-" + subj
    # load the right files from derivatives folder
    derivatives_root = os.path.join(bids_root, "derivatives")
    preproc_folder = os.path.join(derivatives_root, subj, "preprocessing")

    # the files of interest: processing flag = "tsssTransEve"
    bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, root=derivatives_root, processing="tsssTransEve")
    print(f"\n\nThe following files with processing= \"tsssTransEve\" were found: {bids_derivatives.match()}\n\n")
    raw = read_raw_bids(bids_path=bids_derivatives)

# preprocessing
    # filter
    if do_filter:
        raw = prepper.filter_raw(raw, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design, n_jobs=n_jobs)
    
    # resample
    if do_resample:
        raw = prepper.resample_raw(raw, s_freq=s_freq)
    
    # ECG artifacts
    if do_ecg_correction_ssp:
        ecg_artifact = mne.preprocessing.create_ecg_epochs(raw, ch_name=ecg_channel)
        fig, ax = plt.subplots(2,1, figsize=(15,10))
        fig.suptitle("ECG artifacts")
        ecg_artifact.apply_basline((None, None))
        ecg_projs, _ = mne.preprocessing.compute_proj_ecg(raw, n_grad=n_grad, n_mag=n_mag, n_eeg=n_eeg)
        raw.add_proj(ecg_projs, info=raw.info)
        fig.subplot(2,1,1)
        plt.title("ECG-Topomaps")
        ecg_artifact.plot_joint()
        plt.subplot(2,1,2)
        mne.viz.plot_projs_topomap(ecg_projs, info=raw.info)
        savename = os.path.join(preproc_folder, "ECG Topomap")
        plt.savefig(fig, savename)

    #EOG artifacts    
    if do_eog_correction_ssp:
        eog_evoked = mne.preprocessing.create_eog_epochs(raw, eog_channel=eog_channel).average()
        eog_evoked.apply_baseline((None, None))
        eog_projs, _ = mne.preprocessing.compute_proj_eog(raw, n_grad=n_grad, n_mag=n_mag, n_eeg=n_eeg, no_proj=True)
        raw.add_proj(eog_projs, info=raw.info)
        fig = eog_evoked.plot_joint()
        savename = os.path.join(preproc_folder, "EOG Topomap")
        plt.savefig(fig, savename)


########### to do:
# update BIDS with processing-flag + save




# Epochs

    if do_source_loc:
        epochs = mne.Epochs(raw, tmin=-1.5, tmax=1, baseline=(-1.5,-1), on_missing = "ignore")
        print(f"Epochs metadata = {epochs.metadata}")





        epochs_filename = fnr.get_epochs_file(subj, raw_name)
        if event_dict is not None and isinstance(event_dict, dict):
            if not os.path.isfile(epochs_filename):
                raw = mne.io.read_raw_fif(raw_name, preload=True)
                events = mne.find_events(raw)
                epochs = mne.Epochs(raw, events, event_dict, tmin=-1.5, tmax=1, baseline=(-1.5,-1), 
                                                #reject=reject, 
                                                on_missing = 'ignore')
                print(f"Epochs metadata = {epochs.metadata}")
                epochs.save(epochs_filename, overwrite=True)
            else:
                epochs = mne.read_epochs(epochs_filename)
                print(f"Epochs metadata = {epochs.metadata}")
            try:
                del event_dict
            except Exception as e:
                print(f"Something went wrong while epoching {raw_name}")





# concatenate

if concat_raws == True:
    for subj in subjects:
        raw_files = fnr.get_tsss_eve_fifs(subj)
        if len(raw_files) > 1:
            print(f"Concatenating the following files:\n {[raw_files]}")
            concat_name = fnr.get_filename(subj, "concat")
            if not os.path.isfile(concat_name):
                raws = []
                for idx, raw in enumerate(raw_files):
                    rawname = "raw-" + str(idx)
                    rawname = mne.io.read_raw_fif(raw, preload=True)
                    raws.append(rawname)

                raw_concat = mne.concatenate_raws(raws)   # should also work with raw_concat = mne.concatenate_raws([raw_files])???

                if do_filter:
                    raw_concat = prepper.filter_raw(raw_concat, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design, n_jobs=n_jobs)
                if do_resample:
                    raw_concat = prepper.resample_raw(raw_concat, s_freq=s_freq)
                raw_concat.save(concat_name, overwrite=True)   # This should become a BIDS compatible file save later
                print(f"\nConcatenated file info:\n{raw_concat.info}")
            else:
                print("\nNot doing file concatenation, as concat-file exists...")
        else:
            print("Not doing concatenation, as less than 2 fif-files have been found...")
            concat_name = fnr.get_filename(subj, "concat")
            if not os.path.isfile(concat_name):                # Test me!!!!
                try:
                    for raw in raw_files:
                        raw = mne.io.read_raw_fif(raw, preload=True)
                        raw.save(concat_name)
                        break
                except Exception as e:
                    print(f"Attempt to copypaste rawfile to concat-file failed...\n{e}")



"""
To do:

if multiple files, when loading them there should be a run flag --> has to be added

save raws with added events in a BIDS compatible manner
save concat-file in a BIDS compatible manner
ecg_artifact_correction
eog_artifact_correction


"""
