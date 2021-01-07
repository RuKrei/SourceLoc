#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from configuration import (subjects, n_jobs, bids_root, data_root, session, 
                            concat_raws, l_freq, h_freq, fir_design, s_freq,
                            do_filter, do_resample, do_ecg_correction_ssp, 
                            do_ecg_correction_ica, do_ecg_correction_regression,
                            do_eog_correction_ssp, do_eog_correction_ica,
                            do_eog_correction_regression)
import glob
import mne
from utils.utils import FileNameRetriever , RawPreprocessor
import utils.utils as u
import pandas as pd
from shutil import copyfile
from mne_bids import BIDSPath, write_raw_bids, read_raw_bids

fnr = FileNameRetriever(bids_root)
prepper = RawPreprocessor()


# Combine eventfile + raw

for subj in subjects:
    if not subj.startswith("sub-"):
        subj = "sub-" + subj
    subj_root = os.path.join (data_root, subj)
    raws = glob.glob(subj_root + "/data/*_tsss.fif")
    print(f"raw_files: {raws}")
    event_folder = fnr.get_filename(subj, "event_folder")
    for raw in raws:
        event_file = prepper.get_event_file(event_folder, raw)
        print(f"eventfile loaded is {event_file}")
        raw_base_file = str(raw)
        raw_name = str(raw)
        raw_name = str(raw_name.split(".")[0] + "-eve.fif" )   # needs to become the BIDS-name from below later
        raw_name = raw_name.split("/")[-1]
        meg_dir = fnr.get_filename(subj, "meg")
        raw_name = os.path.join(meg_dir, raw_name)
        raw = mne.io.read_raw_fif(raw, preload=True)
        event_dict = None
        if not os.path.isfile(raw_name):
            if os.path.isfile(str(event_file)):
                event_file, event_dict = prepper.transform_eventfile(event_file)
                print(f"raw_file: {raw_name}")
                print(f"event_file: {event_file}")   
                print(raw.info)
                raw.add_events(event_file) 
                print(f"After adding the events --> {raw.info}")
            else:
                print(f"No eventfile found for {raw_base_file}")
            if do_filter:
                raw = prepper.filter_raw(raw, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design, n_jobs=n_jobs)
            if do_resample:
                raw = prepper.resample_raw(raw, s_freq=s_freq)
            
            raw.save(raw_name, overwrite=True)   # This should become a BIDS compatible file save later

            bids_path = BIDSPath(subject=subj, session=session, task="resting", root=bids_root)
            fbase = os.path.join(bids_root, "derivatives", "sub-" + subj)
            fmeg = os.path.join(fbase, "meg")

            newname = bids_path.update(processing="eve")
            newname = newname.basename + ".fif"
        #    write_raw_bids(raw, bids_path, overwrite=True)    # this breaks with fif-files due to an IAS-Key Error
            
        else:
            raw = mne.io.read_raw_fif(raw_name, preload=True)
            if os.path.isfile(str(event_file)):
                event_file, event_dict = prepper.transform_eventfile(event_file)
                #if event_dict == None:
                #    event_dict = "foo"
                print(f"raw_file: {raw_name}")
                print(f"event_file: {event_file}")
            else:
                print(f"No eventfile found for {raw_base_file}")


# Epochs

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

save raws with added events in a BIDS compatible manner
save concat-file in a BIDS compatible manner
ecg_artifact_correction
eog_artifact_correction


"""
