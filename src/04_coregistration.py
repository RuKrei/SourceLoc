# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from os.path import isfile
from configuration import subjects, subjects_dir, derivatives_root, concat_raws, use_single_transfile, session
import mne
from utils.utils import FileNameRetriever
from mne_bids import BIDSPath, read_raw_bids
import glob

fnr = FileNameRetriever(derivatives_root)

for subj in subjects:
    subsubj = "sub-" + subj
    subjects_dir = fnr.get_filename(subsubj, "subjects_dir")
    bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", 
                                root=derivatives_root, processing="tsssTransEvePreproc")
    print(f"\n\nThe following files with processing= \"tsssTransEvePreproc\" were found: {bids_derivatives.match()}\n\n")
    
    if bids_derivatives.match() == []:
        bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", 
                                root=derivatives_root, processing="tsssTransNoEvePreproc")
    
    #load data
    # This should actually be:
    #try:
    #    raw = read_raw_bids(bids_path=bids_derivatives)
    #except Exception as e:
    #    print(e)
    #    print("Looking for a proprocessed fif-file without Events...")
    #    bids_derivatives.update(processing="tsssTransNoEvePreproc")
    
    # loading manually, as BIDS query returns all kinds of things...
    target_dir = os.path.join(derivatives_root, subsubj, "ses-resting", "meg", subsubj)
    try:
        rawfile = glob.glob(target_dir + "*tsssTransEvePreproc_meg.fif")[0]              #### breaks if rawfile is not found!
    except Exception as e:
        print(e)
    # if no events
    if not os.path.isfile(rawfile):
        rawfile = glob.glob(target_dir + "*tsssTransNoEvePreproc_meg.fif")[0]

    
    if use_single_transfile == True:
        transfile = fnr.get_single_trans_file(subsubj)
        if isfile(transfile):
            print(f"Skipping coregistration, because a transfile ({transfile}) already exists")
        else:
            print(f"\n\n\n--> Transfile should be called: {transfile}\n\n\n")
            mne.gui.coregistration(subject=subsubj, subjects_dir=subjects_dir, inst=rawfile) # BIDS: inst=raw.filenames[0])


"""
    
    if concat_raws:
        concat_file = fnr.get_filename(subj=subj, file="concat")
        transfile = fnr.get_trans_file(subj=subj, fif=concat_file)
        subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
        if os.path.isfile(transfile):
            print("\n\n\nCoregistration skipped, as concat-transfile exists.\n\n\n")
        else:
            
    
    fifs = fnr.get_tsss_eve_fifs(subj)
    for fif in fifs:
        if use_single_transfile:
            transfile = fnr.get_single_trans_file(subj=subj)
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


# To do:
# - if "tsssTrans" in processing --> one transfile
# - if "raw" in processing --> specific transfile
#
# maybe only use processing="tsssTransEvePreproc" for frequency spectrum at first
#

"""