# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import numpy as np
from configuration import (subjects, n_jobs, bids_root, use_source_model_for_freq, 
                            pick_meg, pick_eeg, freq_bands, concat_raws)
import mne
from utils.utils import FileNameRetriever, RawPreprocessor

fnr =FileNameRetriever(bids_root)
prepper = RawPreprocessor()

for subj in subjects:
    rawfile = fnr.get_filename(subj, "concat")
    raw = mne.io.read_raw(rawfile, preload=True)
    filebase = concat_file.split("/")[-1].split(".")[0]
    subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
    src = fnr.get_filename(subj=subj, file=use_source_model_for_freq)
    bem_sol = fnr.get_filename(subj=subj, file="3-layer-BEM-sol")
    trans = fnr.get_trans_file(subj, concat_file)
    fwd_name = fnr.get_filename(subj=subj, file="fwd")
    event_folder = fnr.get_filename(subj, "event_folder")
    try:
        event_file = prepper.get_event_file(event_folder, raw)
        event_file, event_dict = prepper.transform_eventfile(event_file)
        raw = prepper.combine_events_w_raw(raw, event_file)
    except Exception as e:
        print(f"\nSomething went wrong for: \n{raw} when trying to locate the Event-File")

    if os.path.isfile(fwd_name):
        fwd = mne.read_forward_solution(fwd_name)
    else:    
        fwd = mne.make_forward_solution(raw.info, src=src, bem=bem_sol,
                                trans=trans, 
                                meg=pick_meg, eeg=pick_eeg, mindist=0.2, 
                                ignore_ref=False, 
                                n_jobs=n_jobs, verbose=True)
        mne.write_forward_solution(fwd_name, fwd)
    
    eve, event_dict = prepper.transform_event_file(eve_file, results_dir)



"""

new_eve_filename = os.path.join(results_dir, eve_file_base + '_processed.csv')
events = u.read_pandas_events_csv(new_eve_filename)
event_id = dict()
event_id = u.generate_spike_dict(events)
print(f"\n\n\nEvents as loaded via pandas: {events}")
events = u.drop_names_from_event_df(events)
raw.add_events(events) #, stim_channel='STI014')
raw = u.do_crop_and_filter(raw, results_dir, crop_tmax=None, n_jobs=n_jobs, filebase=filebase)
raw = u.do_resample(raw, down_sfreq)
raw, fig = u.do_ecg_reduction(raw, filename=filename, for_report_dir=for_report_dir)
fname = (filebase + '_active_ssp_projections_ecg.png')
fname = os.path.join(for_report_dir, fname)
fig.savefig(fname)
#raw = u.do_eog_reduction(raw, show_plot=show_plot, filename=filename, for_report_dir=for_report_dir)                      
raw.save(raw_name, overwrite=True)

"""


"""
To do:

- eLORETA
- dSPM Vector solution
- ECD
"""