# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from configuration import (subjects, n_jobs, bids_root, 
                            concat_raws, l_freq, h_freq, fir_design, s_freq,
                            do_filter, do_resample)
import glob
import mne
from utils.utils import FileNameRetriever , RawPreprocessor
import utils.utils as u
import pandas as pd

fnr = FileNameRetriever(bids_root)
prepper = RawPreprocessor()


# Combine eventfile + raw

for subj in subjects:
    raws = fnr.get_tsss_fifs(subj=subj)
    event_folder = fnr.get_filename(subj, "event_folder")
    for raw in raws:
        if raw.split("/")[-1] == "sub-BF28011991_ses-resting_task-resting_run-00_proc-tsss.fif":
            event_file = prepper.get_event_file(event_folder, raw)
            raw_name = str(raw)
            raw_name = str(raw_name.split(".")[0] + "-eve.fif" )
            event_file, event_dict = prepper.transform_eventfile(event_file)
            print(f"event_file: {event_file}")
            raw = mne.io.read_raw_fif(raw, preload=True)   
            print(raw.info)
            raw.add_events(event_file) 
            print(f"After adding the events --> {raw.info}")
            if do_filter:
                raw = prepper.filter_raw(raw, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design, n_jobs=n_jobs)
            if do_resample:
                raw = prepper.resample_raw(raw, s_freq=s_freq)
            raw.save(raw_name, overwrite=True)   # This should become a BIDS compatible file save later




# concatenate

if concat_raws == True:
    for subj in subjects:
        raw_files = fnr.get_tsss_eve_fifs(subj)
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



"""
To do:

save raws with added events in a BIDS compatible manner
save concat-file in a BIDS compatible manner
extract epochs here, average and save accordingly


"""
