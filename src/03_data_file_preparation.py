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
    bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", 
                                root=derivatives_root, processing="tsssTransEve", suffix="meg")
    print(f"\n\nThe following files with processing= \"tsssTransEve\" were found: {bids_derivatives.match()}\n\n")
    all_raws = bids_derivatives.match()
    
    # if no events have been found, the tsssTransEve file will not exist, so:
    if all_raws == []:
        bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", 
                                root=derivatives_root, processing="tsssTrans", suffix="meg")
        all_raws = bids_derivatives.match()
    
    for run, raw in enumerate(all_raws):
        raw = read_raw_bids(raw)
        run += 1

    # preprocessing
        # filter
        if do_filter:
            raw = prepper.filter_raw(raw, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design, n_jobs=n_jobs)

        # resample
        if do_resample:
            raw = prepper.resample_raw(raw, s_freq=s_freq)

        # ECG artifacts
        if do_ecg_correction_ssp:
            # It's smarter to supervise this step (--> look at the topomaps!)
            try:
                ecg_artifact = mne.preprocessing.create_ecg_epochs(raw, ch_name=ecg_channel)
                ecg_projs, _ = mne.preprocessing.compute_proj_ecg(raw, n_grad=n_grad, n_mag=n_mag, 
                                                                  n_eeg=n_eeg, reject=None)
                # lets not do this now......
                raw.add_proj(ecg_projs, remove_existing=False)
                #raw.apply_proj(ecg_projs, verbose=None)
                fig = mne.viz.plot_projs_topomap(ecg_projs, info=raw.info, show=False)
                savename = os.path.join(preproc_folder, "ECG_projs_Topomap.png")
                fig.savefig(savename)
            except Exception as e:
                print(e)
                print("ECG - Atrifact correction failed!")

        #EOG artifacts    
        if do_eog_correction_ssp:
            # It's a bad idea to do this in an automated step
            try:
                eog_evoked = mne.preprocessing.create_eog_epochs(raw).average()
                #eog_evoked.apply_baseline((None, None))
                eog_projs, _ = mne.preprocessing.compute_proj_eog(raw, n_grad=n_grad, n_mag=n_mag, n_eeg=n_eeg, 
                                                            n_jobs=n_jobs, reject=None)
                raw.add_proj(eog_projs, remove_existing=False).apply_proj()
                figs = eog_evoked.plot_joint(show=False)
                for idx, fig in enumerate(figs):
                    savename = os.path.join(preproc_folder, "EOG Topomap_" + str(idx) + ".png")
                    fig.savefig(savename)
            except Exception as e:
                print(e)
                print("EOG - Atrifact correction failed!")

        # save - if events have been found
        events, event_ids = mne.events_from_annotations(raw)
        if len(events) > 0:
            raw_temp = os.path.join(preproc_folder, "temp.fif")
            raw.save(raw_temp, overwrite=True)
            raw = mne.io.read_raw(raw_temp, preload=False)
            bids_derivatives.update(processing="tsssTransEvePreproc", run=run)                       
            write_raw_bids(raw, bids_derivatives, overwrite=True)
        else:
        # save - if no events
            raw_temp = os.path.join(preproc_folder, "temp.fif")
            raw.save(raw_temp, overwrite=True)
            bids_derivatives.update(processing="tsssTransNoEvePreproc", run=run)   
            raw = mne.io.read_raw(raw_temp, preload=False)                    
            write_raw_bids(raw, bids_derivatives, overwrite=True)
        