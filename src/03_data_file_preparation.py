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

fnr = FileNameRetriever(bids_root)
prepper = RawPreprocessor()

if concat_raws == True:
    for subj in subjects:
        raw_files = fnr.get_tsss_fifs(subj)
        print([raw_files])
        concat_name = fnr.get_filename(subj, "concat")
        if not os.path.isfile(concat_name):
            raws = []
            for idx, raw in enumerate(raw_files):
                rawname = "raw-" + str(idx)
                rawname = mne.io.read_raw(raw, preload=True)
                raws.append(rawname)

            raw_concat = mne.concatenate_raws(raws)

            if do_filter:
                raw_concat = prepper.filter_raw(raw_concat, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design, n_jobs=n_jobs)
            if do_resample:
                raw_concat = prepper.resample_raw(raw_concat, s_freq=s_freq)
            raw_concat.save(concat_name, overwrite=True)
            print(f"\nConcatenated file info:\n{raw_concat.info}")
        else:
            print("\nNot doing file concatenation, as concat-file exists...")


# Combine eventfile + raw

for subj in subjects:
    raws = fnr.get_concat_fif(subj)
    for raw in raws:
        event_folder = fnr.get_filename(subj, "event_folder")
#        try:
        event_file = prepper.get_event_file(event_folder, raw)
        event_file, event_dict = prepper.transform_eventfile(event_file)    
        raw = prepper.combine_events_w_raw(raw, event_file)      ########breaks here!
#        except Exception as e:
#            print(f"\nSomething went wrong for: \n{raw}")
