# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import glob
import mne
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import shutil
import subprocess


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
            Report-Folder:
                "report"
            
            

            
        """
        
        self.subj = subj
        self.file = file
        

        # basic folder structure
        fbase = os.path.join(self.bids_root, subj)
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
            tsss_dir = os.path.join(self.bids_root, + subj, "meg")
            concat_fname = subj + "-concat-raw-tsss.fif"
            concat = os.path.join(tsss_dir, concat_fname)
            return concat

         # filepath - Subjects Directory
        if file == "subjects_dir":
            return os.path.join(self.bids_root, subj, "freesurfer")
        
        # filepath - Frequency distribution - MNE
        if file=="freqMNE":
            return os.path.join(freq, "MNE")
        
        # filepath - Forward solution
        if file=="fwd":
            fwdname = str(subj + "-forward.fif")
            return os.path.join(ffwd, fwdname)
        
        # filepath - Spike Folder
        if file=="spikes":
            fbase = os.path.join(self.bids_root, subj)
            spikes = os.path.join(fbase, "spikes")
            return spikes

        # filepath - Report Folder
        if file == "report":
            return os.path.join(fbase, "report")

        # filepath - Preprocessing Folder
        if file == "preprocessing":
            return os.path.join(fbase, "preprocessing")


    def get_tsss_fifs(self, subj=None):
        tsss_dir = os.path.join(self.bids_root, subj, "meg")
        fifs = glob.glob(tsss_dir + "/*tsss.fif")
        return fifs
    
    def get_tsss_eve_fifs(self, subj=None):
        tsss_dir = os.path.join(self.bids_root, subj, "meg")
        fifs = glob.glob(tsss_dir + "/*tsss-eve.fif")
        return fifs
    
    def get_concat_fif(self, subj=None):
        tsss_dir = os.path.join(self.bids_root, subj, "meg")
        fifs = glob.glob(tsss_dir + "/*concat-raw-tsss.fif")
        return fifs
    
    def get_trans_file(self, subj=None, fif=None):
        trans_dir = os.path.join(self.bids_root, subj, "meg")
        trans_name = fif.split("/")[-1].split(".")[0] + "-trans.fif"
        trans_name = os.path.join(trans_dir, trans_name)
        return trans_name
    
    def get_single_trans_file(self, subj=None):
        trans_dir = os.path.join(self.bids_root, subj, "trans_files")
        trans_name = str(subj) + "-trans.fif"
        trans_name = os.path.join(trans_dir, trans_name)
        return trans_name

    def get_epochs_file(self, subj=None, fif=None): # becomes processing=epo
        epo_dir = os.path.join(self.bids_root, subj, "meg")
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
    
    def filter_raw(self, raw, l_freq, h_freq, n_jobs=1, fir_design="firwin"):
        raw.load_data()
        raw = raw.filter(l_freq=l_freq, h_freq=h_freq, 
                        n_jobs=n_jobs, fir_design=fir_design)
        print("Filtering complete!")
        return raw
    
    def resample_raw(self, raw, s_freq=300, events=None, n_jobs=1):
        print(f"Resampling to {s_freq} Hz")
        raw = raw.resample(s_freq, npad='auto', events=events, n_jobs=n_jobs)
        print("Resampling complete!")
        return raw

    def combine_raw_and_eve(self, rawfile=None, run=1):
        eve_name = rawfile.split(".fif")[0] + "_Events.csv"
        if not os.path.isfile(eve_name):
            eve_name = rawfile.split(".fif")[0] + "_Events.txt"
        if os.path.isfile(eve_name): # if fif-file matches event-file --> add events to fif-file
            print(f"\n\nNow adding Events ({eve_name}) to fif ({rawfile})\n\n")
            raw = mne.io.read_raw(rawfile, preload=True, on_split_missing="ignore")
            event_file, _ = self.transform_eventfile(eve_name)
            raw.add_events(event_file)
            return raw
    
    def raw_to_epoch(self, rawfile=None):
        eve_name = rawfile.split(".fif")[0] + "_Events.csv"
        if not os.path.isfile(eve_name):
            eve_name = rawfile.split(".fif")[0] + "_Events.txt"
        if os.path.isfile(eve_name): # if fif-file matches event-file --> add events to fif-file
            try:
                print(f"\n\nNow epoching events from {rawfile}\n\n")
                event_file, event_dict = self.transform_eventfile(eve_name)
                print(f"\n\nevent_file: {event_file}")
                print(f"\n\nevent_dict: {event_dict}")
                raw = mne.io.read_raw(rawfile)
                epochs = mne.Epochs(raw, events=event_file,
                                        event_id=event_dict, 
                                        tmin=-1.5, tmax=1, 
                                        baseline=(-1.5,-1), on_missing = "ignore", 
                                        event_repeated="merge")
                del(raw)
                return epochs
            except Exception as e:
                print(f"failed at returning an epochs object for: {rawfile}\nbecause of {e}")













def plot_freq_band_dors(stc_band, band=None, subject=None, subjects_dir=None, 
                        filebase=None):
    title = str(filebase) + " " + str(band)
    brain = stc_band.plot(subject=subject, subjects_dir=subjects_dir, 
                        hemi='both',
                        time_label=title, colormap='inferno', 
                        add_data_kwargs=dict(time_label_size=10),
                        clim=dict(kind='percent', lims=(25, 70, 99)))
    brain.show_view(azimuth=0, elevation=0, roll=0)
    return brain

def plot_freq_band_lat(stc_band, band=None, subject=None, subjects_dir=None, filebase=None):
    title = str(filebase) + " " + str(band)
    brain_lh = stc_band.plot(subject=subject, subjects_dir=subjects_dir, hemi='lh',
                        time_label=title, colormap='inferno', size=(1500, 800),
                        add_data_kwargs=dict(time_label_size=10),
                        clim=dict(kind='percent', lims=(25, 70, 99)))
    brain_rh = stc_band.plot(subject=subject, subjects_dir=subjects_dir, hemi='rh',
                        time_label=title, colormap='inferno', size=(1500, 800),
                        add_data_kwargs=dict(time_label_size=10),
                        clim=dict(kind='percent', lims=(25, 70, 99)))                    
    return (brain_lh, brain_rh)

def plot_freq_band_med(stc_band, band=None, subject=None, subjects_dir=None, filebase=None):
    title = str(filebase) + " " + str(band)
    brain_lh = stc_band.plot(subject=subject, subjects_dir=subjects_dir, hemi='lh',
                        views='medial',
                        time_label=title, colormap='inferno', size=(1500, 800),
                        add_data_kwargs=dict(time_label_size=10),
                        clim=dict(kind='percent', lims=(25, 70, 99)))
    brain_rh = stc_band.plot(subject=subject, subjects_dir=subjects_dir, hemi='rh',
                        time_label=title, colormap='inferno', size=(1500, 800),
                        add_data_kwargs=dict(time_label_size=10),
                        views='medial',
                        clim=dict(kind='percent', lims=(25, 70, 99)))                    
    return (brain_lh, brain_rh)
 
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

def run_shell_command(command):
    subprocess.run(command, shell=True, 
                   capture_output=True, check=True)
