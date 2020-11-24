# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import numpy as np
from configuration import (subjects, n_jobs, bids_root, use_source_model_for_sourceloc, 
                            pick_meg, pick_eeg, concat_raws)
import mne
from utils.utils import FileNameRetriever, RawPreprocessor

fnr =FileNameRetriever(bids_root)
prepper = RawPreprocessor()

for subj in subjects:
    concat_file = fnr.get_filename(subj, "concat")
    raw = mne.io.read_raw_fif(concat_file, preload=True)
    filebase = concat_file.split("/")[-1].split(".")[0]
    subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
    src = fnr.get_filename(subj=subj, file=use_source_model_for_sourceloc)
    bem_sol = fnr.get_filename(subj=subj, file="3-layer-BEM-sol")
    trans = fnr.get_trans_file(subj, concat_file)
    fwd_name = fnr.get_filename(subj=subj, file="fwd")
    event_folder = fnr.get_filename(subj, "event_folder")
    event_file = fnr.get_event_file(subj, concat_file)
    print(f"eventfile should be named: {event_file}")
    try:
        event_file, event_dict = prepper.transform_eventfile(event_file)
    except Exception as e:
        print(f"\nSomething went wrong for: \n{raw} when trying to locate the Event-File")

    epochs_filename = fnr.get_epochs_file(subj, concat_file)
    if not os.path.isfile(epochs_filename):
        events = mne.find_events(raw)
        epochs = mne.Epochs(raw, events, event_dict, tmin=-2, tmax=2, baseline=(-2,-1), 
                                        #reject=reject, 
                                        on_missing = 'ignore')
        epochs.save(epochs_filename, overwrite=True)
    else:
        epochs = mne.read_epochs(epochs_filename)

    if os.path.isfile(fwd_name):
        fwd = mne.read_forward_solution(fwd_name)
    else:    
        fwd = mne.make_forward_solution(epochs.info, src=src, bem=bem_sol,
                                trans=trans, 
                                meg=pick_meg, eeg=pick_eeg, mindist=0.2, 
                                ignore_ref=False, 
                                n_jobs=n_jobs, verbose=True)
        mne.write_forward_solution(fwd_name, fwd)

    noise_cov = mne.compute_covariance(epochs, tmax=-1, 
                                        #projs=, 
                                        method='auto',
                                        n_jobs=n_jobs)
    data_cov = mne.compute_covariance(epochs,
                                        tmin=-0.5, 
                                        tmax=0.3, 
                                        #projs=, 
                                        method='auto',
                                        n_jobs=n_jobs)

    eventnames = epochs.event_id.keys()
    for eventname in eventnames:
        if str(eventname) == 'ignore_me' or str(eventname) == 'AAA':
            pass
        else:
            save_dir = fnr.get_filename(subj, "spikes")
            save_name = str(subj + "_" + eventname + '-ave.fif')
            save_name = os.path.join(save_dir, eventname, save_name)
            if not os.path.isfile(save_name):
                plot_event = epochs[str(eventname)].load_data().average()
                current_event_dir = save_name.split("/")[0:-2]
                if not os.path.isdir(current_event_dir):
                    os.mkdir(current_event_dir)
                    os.mkdir(os.path.join(current_event_dir, 'custom_pics'))
                    os.mkdir(os.path.join(current_event_dir, 'custom_time_series'))
                plot_event.save(save_name)
            else:
                plot_event = mne.read_epochs(save_name)


"""
To do:

- eLORETA
- dSPM Vector solution
- ECD
"""