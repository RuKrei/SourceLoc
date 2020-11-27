# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import glob
import mne
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder


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
            Concatednated -trans-tsss.fif-file:
                "concat"
            Eventfolder:
                "event_folder"
            Freesurfer Subjects Directory:
                "subjects_dir"
            Frequency distribution - MNE:
                "freqMNE"
            Forward solution:
                "fwd"
            Spikes:
                "spikes"
            MEG-Folder:
                "meg"
            
            

            
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
        
        # filepath - concat file
        if file == "concat":
            tsss_dir = os.path.join(self.bids_root, "derivatives", "sub-" + subj, "meg")
            concat_fname = subj + "-concat-raw-tsss.fif"
            concat = os.path.join(tsss_dir, concat_fname)
            return concat
        
        # filepath - event file
        if file == "event_folder":
            return os.path.join(fbase, "eventfiles")

         # filepath - Subjects Directory
        if file == "subjects_dir":
            return os.path.join(self.bids_root, "derivatives", "sub-" + subj, "freesurfer")
        
        # filepath - Frequency distribution - MNE
        if file=="freqMNE":
            return os.path.join(freq, "MNE")
        
        # filepath - Frequency distribution - DICS
        if file=="freqDICS":
            return os.path.join(freq, "DICS")
        
        # filepath - Forward solution
        if file=="fwd":
            fwdname = str(subj + "-forward.fif")
            return os.path.join(ffwd, fwdname)
        
        # filepath - Spike Folder
        if file=="spikes":
            fbase = os.path.join(self.bids_root, "derivatives", "sub-" + self.subj)
            spikes = os.path.join(fbase, "spikes")
            return spikes

        # filepath - meg
        if file=="meg":
            fbase = os.path.join(self.bids_root, "derivatives", "sub-" + self.subj)
            megfolder = os.path.join(fbase, "meg")
            return megfolder


    def get_tsss_fifs(self, subj=None):
        tsss_dir = os.path.join(self.bids_root, "derivatives", "sub-" + subj, "meg")
        fifs = glob.glob(tsss_dir + "/*tsss.fif")
        return fifs
    
    def get_tsss_eve_fifs(self, subj=None):
        tsss_dir = os.path.join(self.bids_root, "derivatives", "sub-" + subj, "meg")
        fifs = glob.glob(tsss_dir + "/*tsss-eve.fif")
        return fifs
    
    def get_concat_fif(self, subj=None):
        tsss_dir = os.path.join(self.bids_root, "derivatives", "sub-" + subj, "meg")
        fifs = glob.glob(tsss_dir + "/*concat-raw-tsss.fif")
        return fifs
    
    def get_trans_file(self, subj=None, fif=None):
        trans_dir = os.path.join(self.bids_root, "derivatives", "sub-" + subj, "meg")
        trans_name = fif.split("/")[-1].split(".")[0] + "-transfile.fif"
        trans_name = os.path.join(trans_dir, trans_name)
        return trans_name
    
    def get_single_trans_file(self, subj=None):
        trans_dir = os.path.join(self.bids_root, "derivatives", "sub-" + subj, "meg")
        trans_name = str(subj) + "-transfile.fif"
        trans_name = os.path.join(trans_dir, trans_name)
        return trans_name

    def get_event_file(self, subj=None, fif=None):
        fbase = os.path.join(self.bids_root, "derivatives", "sub-" + subj)
        event_dir = os.path.join(fbase, "eventfiles")
        eve_name = fif.split("/")[-1].split(".")[0] + "-eve.csv"
        if not os.path.isfile(eve_name):
            eve_name = fif.split("/")[-1].split(".")[0] + "-eve.txt"
        eve_name = os.path.join(event_dir, eve_name)
        return eve_name

    def get_epochs_file(self, subj=None, fif=None):
        epo_dir = os.path.join(self.bids_root, "derivatives", "sub-" + subj, "meg")
        epo_name = fif.split("/")[-1].split(".")[0] + "-epo.fif"
        epo_name = os.path.join(epo_dir, epo_name)
        return epo_name

class RawPreprocessor():

    def __init__(self):
        pass

    def transform_eventfile(self, event_file):
        """
        Receives a .csv or .txt files as exported i.e. via brainstorm and
        transforms it according to mne-needs.

        Returns a tuple:
        new_eve_file, event_dict
        --> = transformed event-file, dictionary with Eventnames
        """

        eve = pd.read_csv(event_file, header=0)
        le = LabelEncoder()
        labels = eve.iloc[:,0]
        print(f"Labels --> {labels}")
        l_enc = le.fit_transform(labels)
        l_enc = l_enc
        new_eve_file = pd.DataFrame([eve.iloc[:,1], eve.iloc[:,0], (l_enc +1)]).T
        new_eve_file.reset_index(drop=True, inplace = True)
        new_eve_file.iloc[:,0] = (new_eve_file.iloc[:,0]*1000).astype(int)
        new_eve_file.iloc[0,2] = 0  #create one pseudo-event (that is going to be dropped later for some reason)
        new_eve_file.iloc[:,0] = new_eve_file.iloc[:,0].astype(int)
        new_eve_file.iloc[:,1] = new_eve_file.iloc[:,1].astype(str)
        new_eve_file.iloc[:,2] = new_eve_file.iloc[:,2].astype(int)
        new_eve_file.iloc[:,1] = int("0")

        name_of_events = np.unique(eve.iloc[:,0])
        name_of_events = np.sort(name_of_events)
        event_dict=dict()
        event_dict['ignore_me']=0
        for i in range(name_of_events.size):
            key = (name_of_events[i])
            val = i + 1
            event_dict[key] = val
        return new_eve_file, event_dict
    
    def combine_events_w_raw(self, raw, event_file):
        """
        Receives the eventfile as processed by RawPreprocessor.transform_eventfile() + loaded raw file.
        Returns raw with events added
        """
        print(f"event_file --> {event_file.iloc[:,0]}")
        raw = raw.add_events(event_file)
        return raw

    def get_event_file(self, event_folder, raw_filename):
        rawfile = raw_filename.split("/")[-1].split(".")[0]
        eventfile = os.path.join(event_folder, rawfile + "_Events.csv")
        if os.path.isfile(eventfile):
            return eventfile
        eventfile = os.path.join(event_folder, rawfile + "_Events.txt")
        if os.path.isfile(eventfile):
            return eventfile
    
    def filter_raw(self, raw, l_freq=1, h_freq=70, fir_design="firwin", n_jobs=1):
        raw = raw.filter(l_freq=l_freq, h_freq=h_freq, fir_design=fir_design)
        return raw
    
    def resample_raw(self, raw, s_freq=300):
        raw.resample(s_freq, npad='auto')
        return raw

def plot_freq_band_dors(stc_band, band=None, subject=None, subjects_dir=None, filebase=None):
    title = (filebase + ' - Frequenzanalyse - ' + band)
    brain = stc_band.plot(subject=subject, subjects_dir=subjects_dir, hemi='both',
                        time_label=title, colormap='inferno', 
                        clim=dict(kind='percent', lims=(25, 70, 99)))
    brain.show_view(dict(azimuth=0, elevation=0), roll=0)
    return brain
def plot_freq_band_lat(stc_band, band=None, subject=None, subjects_dir=None, filebase=None):
    title = (filebase + ' - Frequenzanalyse - ' + band)
    brain = stc_band.plot(subject=subject, subjects_dir=subjects_dir, hemi='split',
                        time_label=title, colormap='inferno', size=(1500, 800),
                        clim=dict(kind='percent', lims=(25, 70, 99)))
    return brain
 
def get_peak_points(stc, hemi='lh', 
                        tmin=-.02, tmax=0, nr_points=5, 
                        mode='abs'):
    '''
    Returns a list of vertices with peak activation at time points 
    as given by tmin, tmax and nr_points
    '''
    pos = dict()
    for t in np.linspace(tmin, tmax, nr_points):
        pos[t], _ = stc.get_peak(hemi=hemi, tmin=(t -0.005), tmax=(t + 0.005), mode=mode)
    return pos.values()