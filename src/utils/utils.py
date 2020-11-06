# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import glob

class FileNameRetriever():
    def __init__(self, bids_root):
        self.bids_root = bids_root

    def get_filename(self, subj=None, file=None):
        """
        Retrive filenames for a subject
        Parameters:
        subj: str
            Subject
        file: str
            The desired file as saved by the pipeline:
            Source models:
                "ico4" , "ico5", "oct5", "oct6", "vol-src"
            BEM models:
                "3-layer-BEM-model", "single-shell-model"
            BEM - Solution:
                "single-shell-BEM-sol", "3-layer-BEM-sol"
        """
        
        self.subj = subj
        self.file = file
        

        # basic folder structure
        fbase = os.path.join(self.bids_root, "derivatives", "sub-" + self.subj)
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
        freport = os.path.join(fbase, "report")

        # filepaths - Source Models
        if file in ["ico4" , "ico5", "oct5", "oct6", "vol-src"]:
            srcfilename = (subj + '-' + file + '-src.fif')
            srcfilename = os.path.join(fsrc, srcfilename)
            return srcfilename
        
        # filepaths - BEM solutions
        if file in ["single-shell-BEM-sol", "3-layer-BEM-sol"]:
            if file == "single-shell-BEM-sol":
                bem_sol_filename = f"{subj}-single-shell-BEM-solution.fif"
                bem_sol_filename = os.path.join(ffwd, bem_sol_filename)
                return bem_sol_filename
            if file == "3-layer-BEM-sol":
                bem_sol_filename = f"{subj}-3-layer-BEM-solution.fif"
                bem_sol_filename = os.path.join(ffwd, bem_sol_filename)
                return bem_sol_filename

        # filepaths - BEM models
        if file in ["3-layer-BEM-model", "single-shell-model"]:
            bem_f_name = subj + '-head.fif'
            if file == "single-shell-model":
                bem_save_name = os.path.join(ffwd, bem_f_name.split(".")[0] + "-single-shell-bemfile-ico4.fif")
                return bem_save_name
            if file == "3-layer-BEM-model":
                bem_save_name = os.path.join(ffwd, bem_f_name.split(".")[0] + "-3-layer-bemfile-ico4.fif")
                return bem_save_name

        def get_tsss_fifs(self, subj=None):
            tsss_dir = os.path.join(self.bids_root, "derivatives", "sub-", subj + "meg")
            fifs = glob.glob(tsss_dir + "*.fif")
            return fifs
