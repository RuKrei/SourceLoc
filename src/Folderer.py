#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from os.path import join as opj
import shutil
import glob
import platform
import argparse


class DerivativesFoldersCreator:
    def __init__(self, BIDS_root, extras_directory, subject):
        self.BIDS_root = BIDS_root
        self.subject = subject
        if not self.subject.startswith("sub-"):
            self.subject = "sub-" + self.subject
        self.fbase = opj(self.BIDS_root, "derivatives", self.subject)
        self.fsrc = opj(self.fbase, "source_model")
        self.fanat = opj(self.fbase, "freesurfer")
        self.fprep = opj(self.fbase, "preprocessing")
        self.spikes = opj(self.fbase, "spikes")
        self.freq = opj(self.fbase, "frequency_distribution")
        self.conn = opj(self.fbase, "connectivity")
        self.ftrans = opj(self.fbase, "trans_files")
        self.freport = opj(self.fbase, "report")
        self.extras_directory = extras_directory
    
    def _create_subfolders(self):
        for fld in [self.fbase, 
                self.fsrc, 
                self.fanat, 
                self.fprep, 
                self.spikes, 
                self.freq, 
                self.conn, 
                self.ftrans, 
                self.freport]:
            if not os.path.exists(fld):
                os.makedirs(fld, exist_ok=True)
                print(f"--> {fld} created.")

    def _recursive_overwrite(self, src, dest, ignore=None):
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
                    self._recursive_overwrite(opj(src, f), 
                                                opj(dest, f), 
                                                ignore)
        else:
            try:
                shutil.copyfile(src, dest)
            except:
                pass
    
    def _copy_data_from_extras(self):
        splitter = "\\" if platform.system().lower().startswith("win") else "/"
        disclaimer = opj(self.extras_directory, "MEG_disclaimer.png")
        title = opj(self.extras_directory, "MEG_title.png")
        report_files = [disclaimer, title]
        ana_dir = opj(self.extras_directory, "anatomy_templates")
        ana_files = glob.glob(ana_dir) 
        for f in report_files:
            try:
                fi = f.split(splitter)[-1]
                target = os.path.join(self.freport, fi)
                shutil.copyfile(f, target)
            except Exception as e:
                print(e)
            #copy fsaverage + fsaverage_sym to local subjects anatomy folder
            for f in ana_files:
                if not (f.split(splitter)[-1].endswith(".png")):
                    try:
                        self._recursive_overwrite(f, self.fanat)
                        print(f"Copying: {f}")
                    except Exception as e:
                        print(e)    
    
    def make_derivatives_folders(self):
        self._create_subfolders()
        if os.path.isdir(self.extras_directory):
            self._copy_data_from_extras()
        

def main():
    parser = argparse.ArgumentParser()     
    parser.add_argument("--BIDS_root", action="store", type=str, required=True, 
                        help="Specify the BIDS root folder")
    parser.add_argument("--subject", action="store", type=str, required=True, 
                        help="Subject name")
    parser.add_argument("--extras_directory", action="store", type=str, required=False, 
                        help="Specify location of extras directory")

    args = parser.parse_args()
    dfc = DerivativesFoldersCreator(args.BIDS_root, args.extras_directory, args.subject)
    dfc.make_derivatives_folders()


if __name__ == '__main__':
    main()