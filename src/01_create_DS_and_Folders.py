# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from mne_bids import BIDSPath, write_raw_bids, print_dir_tree, write_anat
import os
import utils as u
from configuration import subjects, session
import glob
import mne


# BIDS inputs
bids_root = "/home/idrael/DATA/MEG/BIDS_clinic"
data_root = "/home/idrael/DATA/MEG_playground"


# BIDSify
for subj in subjects:
    data_dir = os.path.join(data_root, subj, "data")
    raws = glob.glob(data_dir + "/*tsss.fif")
    for run, rawfile in enumerate(raws):
        run += 1
        bids_path = BIDSPath(subject=subj, session=session, run=run, task="resting", root=bids_root)
        raw = mne.io.read_raw(rawfile)
        write_raw_bids(raw, bids_path, overwrite=True)


# Add derivatives folder structure
for subj in subjects:
    fbase = os.path.join(bids_root, "derivatives", "sub-" + subj)
    fanat = os. path.join(fbase, "freesurfer")
    ffwd = os.path.join(fbase, "forward_model")
    fcss = os.path.join(ffwd, "cortex_source_space")
    fvss = os.path.join(ffwd, "volume_source_space")
    finv = os.path.join(fbase, "inverse_model")
    fci = os.path.join(finv, "cortex_inverse")
    fvi = os.path.join(finv, "volume_inverse")
    spikes = os.path.join(fbase, "spikes")
    freq = os.path.join(fbase, "frequency_distribution")
    fDICS = os.path.join(freq, "DICS")
    fMNE = os.path.join(freq, "MNE")
    conn = os.path.join(fbase, "connectivity")
    cnonspike = os.path.join(conn, "nonspike_all_to_all")
    cspike = os.path.join(conn, "spike")    # this needs to be filled according to spike name later
    feve = os.path.join(fbase, "eventfiles")
    freport = os.path.join(fbase, "report")
    
    folder_list = [fbase, fanat, ffwd, fcss, fvss, finv, fci, fvi, spikes, freq, fDICS, fMNE, conn, cnonspike, cspike, feve, freport]

    for fld in folder_list:
        if not os.path.exists(fld):
            os.makedirs(fld, exist_ok=True)
  


def find_dcm_in_dir(s):
    do_recon_dict = dict()
    try:
        counter = 1
        for file in glob.glob((s + '/1*/**/1000*'), recursive=True):
            do_recon_dict[s] = file
            counter += 1
            if counter > 4:
                return file
                break
    except OSError as e:
        print(e)










print_dir_tree(bids_root)