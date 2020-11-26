# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import numpy as np
from configuration import (subjects, n_jobs, bids_root, use_source_model_for_sourceloc, 
                            pick_meg, pick_eeg, concat_raws)
import mne
from utils.utils import FileNameRetriever, RawPreprocessor
import glob

fnr =FileNameRetriever(bids_root)
prepper = RawPreprocessor()

for subj in subjects:
    meg_folder = fnr.get_filename(subj, "meg")
    epos = glob.glob(meg_folder + "*.epo.fif")
    for epo in epos:
        epochs = mne.read_epochs(epo, preload=True)
        print(f"Event-IDs are: {epochs.event_id}")
    


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