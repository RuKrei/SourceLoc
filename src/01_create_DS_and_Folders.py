# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from mne_bids import BIDSPath, write_raw_bids, print_dir_tree, write_anat
import os
import utils as u
from configuration import subjects, session
import glob
import mne
from dicom2nifti import convert_directory


# BIDS inputs
bids_root = "/home/idrael/DATA/MEG/BIDS_clinic"
data_root = "/home/idrael/DATA/MEG_playground"


def convert_dcm_folder(subj):
    try:
        anafolder = os.path.join(data_root, subj, "data", "anat", subj)
        print(f"Anafolder: {anafolder}")
        folder = str(glob.glob((anafolder + '/1*/100*/100*'), recursive=True)[0])
        print(f"Folder = {folder}")
        convert_directory(folder, anafolder, compression=True, reorient=True)
    except OSError as e:
        print(e)


# BIDSify
for subj in subjects:
    data_dir = os.path.join(data_root, subj, "data")
    raws = glob.glob(data_dir + "/*tsss.fif")
    bids_path = BIDSPath(subject=subj, session=session, task="resting", root=bids_root)
    try:
        for run, rawfile in enumerate(raws):
            run += 1
            raw = mne.io.read_raw(rawfile)
            write_raw_bids(raw, bids_path, overwrite=True)   # tSSS files are put in meg folder (and NOT in derivatives/... folder)
            # this breaks due to an IAS-Key Error in mne_bids
    except Exception as e:
        print(e)
    convert_dcm_folder(subj)
    anafolder = os.path.join(data_root, subj, "data", "anat", subj)
    niis = glob.glob(anafolder + "/*.nii*")
    try:
        for n in niis:
            write_anat(n, bids_path=bids_path, overwrite=True)
    except Exception as e:
        print(e)




# Add derivatives folder structure
for subj in subjects:
    fbase = os.path.join(bids_root, "derivatives", "sub-" + subj)
    fanat = os. path.join(fbase, "freesurfer")
    fprep = os.path.join(fbase, "preprocessing")
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
    
    folder_list = [fbase, fanat, fprep, ffwd, fcss, fvss, finv, fci, fvi, spikes, freq, fDICS, fMNE, conn, cnonspike, cspike, feve, freport]

    for fld in folder_list:
        if not os.path.exists(fld):
            os.makedirs(fld, exist_ok=True)
  

print_dir_tree(bids_root)