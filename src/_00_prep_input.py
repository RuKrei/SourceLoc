#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

"""
This file prepares the input.
It expects a folder with the subject-name/code containing DICOM files in subfolder 1*/100*/100*/xyz.dcm
(This is how DICOMs are being extracted from IMPAX in our clinic)
It creates a BIDS compatible folder structure and moves/renames files accordingly
DICOM-directories are converted to .nii.gz before moving

"""


import os
import shutil
from configuration import input_folder, data_root, subjects


# helper functions
def recursive_overwrite(src, dest, ignore=None):
    if os.path.isdir(src):
        if not os.path.isdir(dest):
            os.makedirs(dest)
        files = os.listdir(src)
        if ignore is not None:
            ignored = ignore(src, files)
        else:
            ignored = set()
        for f in files:
            if f not in ignored:
                recursive_overwrite(os.path.join(src, f), 
                                    os.path.join(dest, f), 
                                    ignore)
    else:
        shutil.copyfile(src, dest)


for s in subjects:
    # define target directory
    target = os.path.join(data_root, s)
    print(f"Target-directory --> {target}")
    if os.path.isdir(target):
        pass
    else:
        # create target directory
        os.mkdir(target)
        print(f"Directory >> {target} << created.")
        # copy files to target directory
        destination = os.path.join(target, s)
        source = os.path.join(input_folder, s)
        try:
            recursive_overwrite(source, destination)
        except FileNotFoundError:
            print(f"No DICOM files found for Patient: {s}")