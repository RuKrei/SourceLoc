#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import shutil
from configuration import input_folder, data_root
import argparse
from os.path import join as opj
from dicom2nifti import convert_directory
import glob

# configuration
bids_root = "./BIDS"
data_root = "./templates"
input_folder = "C:\\Users\\rudik\\MEG\\playground\\input_folder"
open_mp = 4
n_jobs = 4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject", action="store", type=str, required=False)
    parser.add_argument("--bidsroot", action="store", type=str, required=False)
    parser.add_argument("--inputfolder", action="store", type=str, required=False)
    args = parser.parse_args()  

# define subject
    subject = args.subject    
    if not subject:
        poss = [s for s in os.listdir(input_folder)]
        print(f"No subject specified, maybe you want to choose from those:\n {poss}")
        subject = input()
    
    print(f"Subject = {subject}")

# make sure we have an MRI, convert, if necessary
    anafolder = opj(input_folder, subject)
    nii = glob.glob(opj(input_folder, subject, "*.nii.gz"))
    if os.path.isdir(anafolder) and nii == []:
        folder = str(glob.glob((anafolder + '/1*/100*/100*'), recursive=True)[0])
        try:
            convert_directory(
                dicom_directory=folder, 
                output_folder=anafolder, 
                compression=True, 
                reorient=True)
            print("\nMRI was converted to .nii.gz\n")
        except Exception as e:
            print(f"Something went wrong trying to convert the MRI to nifti: {e}")
    nii = glob.glob(opj(input_folder, subject, "*.nii.gz"))
    if nii == []:
        print("No anatomical data found, did you provide an MRI?")
    else:
        print(f"\nMRI already is in .nii.gz-Format: {nii}\nDoing nothing...\n")




if __name__ == '__main__':
    main()