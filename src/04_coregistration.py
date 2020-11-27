# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from configuration import subjects, bids_root, concat_raws, use_single_transfile
import mne
from utils.utils import FileNameRetriever

fnr = FileNameRetriever(bids_root)

for subj in subjects:
    if concat_raws:
        concat_file = fnr.get_filename(subj=subj, file="concat")
        transfile = fnr.get_trans_file(subj=subj, fif=concat_file)
        subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
        if os.path.isfile(transfile):
            print("\n\n\nCoregistration skipped, as concat-transfile exists.\n\n\n")
        else:
            print(f"\n\n\nTransfile should be called: {transfile}\n\n\n")
            mne.gui.coregistration(subject=subj, subjects_dir=subjects_dir, inst=concat_file)
    
    fifs = fnr.get_tsss_eve_fifs(subj)
    for fif in fifs:
        if use_single_transfile:
            transfile = fnr.get_single_trans_file(subj=subj, fif=fif)
            subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
            if os.path.isfile(transfile):
                print(f"\n\n\nCoregistration skipped, as transfile for {fif} exists.\n\n\n")
            else:
                print(f"\n\n\nTransfile should be called: {transfile}\n\n\n")
                mne.gui.coregistration(subject=subj, subjects_dir=subjects_dir, inst=fif)
        
        else:
            transfile = fnr.get_trans_file(subj=subj, fif=fif)
            subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
            if os.path.isfile(transfile):
                print(f"\n\n\nCoregistration skipped, as transfile for {fif} exists.\n\n\n")
            else:
                print(f"\n\n\nTransfile should be called: {transfile}\n\n\n")
                mne.gui.coregistration(subject=subj, subjects_dir=subjects_dir, inst=fif)
