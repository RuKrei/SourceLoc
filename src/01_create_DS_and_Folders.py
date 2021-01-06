#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from mne_bids import BIDSPath, write_raw_bids, print_dir_tree, write_anat, make_dataset_description
import os
import mne_bids
from configuration import subjects, session, bids_root, data_root
from utils.utils import RawPreprocessor
import glob
import mne
from dicom2nifti import convert_directory


def convert_dcm_folder(subj):
    try:
        anafolder = os.path.join(data_root, subj, "data", "anat", subj)
        folder = str(glob.glob((anafolder + '/1*/100*/100*'), recursive=True)[0])
        convert_directory(folder, anafolder, compression=True, reorient=True)
    except Exception as e:
        print(e)

prepper = RawPreprocessor()

# 1. Add derivatives folder structure + bids_root directory, if it doesn't exist 
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
    fMNE = os.path.join(freq, "MNE")
    conn = os.path.join(fbase, "connectivity")
    feve = os.path.join(fbase, "eventfiles")
    ftrans = os.path.join(fbase, "trans_files")
    freport = os.path.join(fbase, "report")
    
    folder_list = [fbase, fmeg, fsrc, fanat, fprep, ffwd, finv, spikes, freq, fMNE, conn, feve, ftrans, freport]

    for fld in folder_list:
        if not os.path.exists(fld):
            os.makedirs(fld, exist_ok=True)


# 2. BIDSify data
for subj in subjects:
    data_dir = os.path.join(data_root, subj)
    raws = glob.glob(data_dir + "/*.fif")
    bids_path = BIDSPath(subject=subj, session=session, task="resting", root=bids_root, processing=None)
    fbase = os.path.join(bids_root, "derivatives", "sub-" + subj)
    fmeg = os.path.join(fbase, "meg")
    for run, rawfile in enumerate(raws):
        if "tsss" in rawfile:
            # --> search for matching eventfile and combine
            eve_name = rawfile.split(".fif")[0] + "_Events.csv"
            if not os.path.isfile(eve_name):
                eve_name = rawfile.split(".fif")[0] + "_Events.txt"
            if os.path.isfile(eve_name):
                print(f"\n\nNow adding Events ({eve_name}) to fif ({rawfile})\n\n")
                raw = mne.io.read_raw(rawfile, preload=True)
                event_file, event_dict = prepper.transform_eventfile(eve_name)
                raw.add_events(event_file)
                bids_path.update(processing="eve")
                raw.save(rawfile, overwrite=True)
        else:
            print("\n\nNo eventfile found for {rawfile}, processing as raw file.\n\n")
            bids_path.update(processing=None)
            # write raw BIDS, as below
        raw = mne.io.read_raw(rawfile)
        mne_bids.write_raw_bids(raw, bids_path, event_id=event_dict, events_data=event_file.to_numpy(), overwrite=True)
    
convert_dcm_folder(subj)
anafolder = os.path.join(data_root, subj, subj)
niis = glob.glob(anafolder + "/*.nii*")
try:
    for n in niis:
        write_anat(n, bids_path=bids_path, overwrite=True)
except Exception as e:
    print(e)

# 3. DataSet description
make_dataset_description(bids_root, name="CDK Epilepsy Dataset", data_license="closed", authors="Rudi Kreidenhuber", overwrite=True)

print_dir_tree(bids_root)



"""

To do:

add disclaimer + title file to report dir --> is accessed from the .extras - Directory

"""