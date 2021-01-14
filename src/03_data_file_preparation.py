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
    subsubj = "sub-" + subj
    # load the right files from derivatives folder
    derivatives_root = os.path.join(bids_root, "derivatives")
    preproc_folder = os.path.join(derivatives_root, subsubj, "preprocessing")

    # the files of interest: processing flag = "tsssTransEve"
    bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", root=derivatives_root, processing="tsssTransEve")
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
        ecg_projs, _ = mne.preprocessing.compute_proj_ecg(raw, n_grad=n_grad, n_mag=n_mag, n_eeg=n_eeg)
        raw.add_proj(ecg_projs)
        fig = mne.viz.plot_projs_topomap(ecg_projs, info=raw.info)
        savename = os.path.join(preproc_folder, "ECG_projs_Topomap.png")
        fig.savefig(savename)

    #EOG artifacts    
    if do_eog_correction_ssp:
        eog_evoked = mne.preprocessing.create_eog_epochs(raw).average()
        #eog_evoked.apply_baseline((None, None))
        eog_projs, _ = mne.preprocessing.compute_proj_eog(raw, n_grad=n_grad, n_mag=n_mag, n_eeg=n_eeg, no_proj=True)
        raw.add_proj(eog_projs)
        figs = eog_evoked.plot_joint()
        for idx, fig in enumerate(figs):
            savename = os.path.join(preproc_folder, "EOG Topomap_" + str(idx) + ".png")
            fig.savefig(savename)

    # save
    events, event_ids = mne.events_from_annotations(raw)
    raw_temp = os.path.join(preproc_folder, "temp.fif")
    raw.save(raw_temp, overwrite=True)
    raw = mne.io.read_raw(raw_temp, preload=False)
    bids_derivatives.update(processing="tsssTransEvePreproc")                       
    write_raw_bids(raw, bids_derivatives, overwrite=True)                           ########### to do: does this include event data/ annotations???

"""
# Epochs

    if do_source_loc:
        epochs = mne.Epochs(raw, events=events, event_id=event_ids, tmin=-1.5, tmax=1, baseline=(-1.5,-1), on_missing = "ignore")
        print(f"Epochs metadata = {epochs.metadata}")
        bids_derivatives.update(processing="tsssTransEpoPreproc")
        write_raw_bids(raw, bids_derivatives, events_data=events, event_id=event_ids, overwrite=True)     



# check if everything worked as expected
    raw = read_raw_bids(bids_derivatives)   
    print(raw.annotations)
"""



"""
To do:
    if multiple files:
    when loading them there should be a run flag --> has to be added
    save a concat-file in a BIDS compatible manner
"""
