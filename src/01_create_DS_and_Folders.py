# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from mne_bids import BIDSPath, write_raw_bids, print_dir_tree, write_anat
import os
from configuration import subjects, session, bids_root, data_root
import glob
from dicom2nifti import convert_directory
from shutil import copyfile


def convert_dcm_folder(subj):
    try:
        anafolder = os.path.join(data_root, subj, "data", "anat", subj)
        folder = str(glob.glob((anafolder + '/1*/100*/100*'), recursive=True)[0])
        convert_directory(folder, anafolder, compression=True, reorient=True)
    except Exception as e:
        print(e)


# Add derivatives folder structure + bids_root directory 
for subj in subjects:
    fbase = os.path.join(bids_root, "derivatives", "sub-" + subj)
    fmeg = os.path.join(fbase, "meg")
    fsrc = os.path.join(fbase, "source_model")
    fanat = os. path.join(fbase, "freesurfer")
    fprep = os.path.join(fbase, "preprocessing")
    ffwd = os.path.join(fbase, "forward_model")
    finv = os.path.join(fbase, "inverse_model")
    spikes = os.path.join(fbase, "spikes")
    freq = os.path.join(fbase, "frequency_distribution")
    fDICS = os.path.join(freq, "DICS")
    fMNE = os.path.join(freq, "MNE")
    conn = os.path.join(fbase, "connectivity")
    cnonspike = os.path.join(conn, "nonspike_all_to_all")
    cspike = os.path.join(conn, "spike")    # this needs to be filled according to spike name later
    feve = os.path.join(fbase, "eventfiles")
    ftrans = os.path.join(fbase, "trans_files")
    freport = os.path.join(fbase, "report")
    
    folder_list = [fbase, fmeg, fsrc, fanat, fprep, ffwd, finv, spikes, freq, fDICS, fMNE, conn, cnonspike, cspike, feve, ftrans, freport]

    for fld in folder_list:
        if not os.path.exists(fld):
            os.makedirs(fld, exist_ok=True)


# BIDSify
for subj in subjects:
    data_dir = os.path.join(data_root, subj, "data")
    raws = glob.glob(data_dir + "/*.fif")
    bids_path = BIDSPath(subject=subj, session=session, task="resting", root=bids_root)
    fbase = os.path.join(bids_root, "derivatives", "sub-" + subj)
    fmeg = os.path.join(fbase, "meg")
    for run, rawfile in enumerate(raws):
        if "sss" in rawfile:
            newname = bids_path.update(processing="tsss", run=run)
            newname = newname.basename + ".fif"
            destination = os.path.join(fmeg, newname)
            copyfile(rawfile, destination)
        #else:
        #    raw = mne.io.read_raw(rawfile)
        #    bids_path.update(run=run)
        #    write_raw_bids(raw, bids_path, overwrite=True)    # this breaks with fif-files due to an IAS-Key Error
    convert_dcm_folder(subj)
    anafolder = os.path.join(data_root, subj, "data", "anat", subj)
    niis = glob.glob(anafolder + "/*.nii*")
    try:
        for n in niis:
            write_anat(n, bids_path=bids_path, overwrite=True)
    except Exception as e:
        print(e)

print_dir_tree(bids_root)