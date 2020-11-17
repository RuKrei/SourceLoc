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
    try:
        event_file, event_dict = prepper.transform_eventfile(event_file)
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





"""
To do:

- eLORETA
- dSPM Vector solution
- ECD
"""