# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os

#from numpy.core.numeric import load
import numpy as np
from configuration import (subjects, n_jobs, bids_root, use_source_model_for_freq, 
                            pick_meg, pick_eeg, freq_bands, concat_raws)
import mne
from mne.minimum_norm import compute_source_psd, make_inverse_operator
from utils.utils import FileNameRetriever, plot_freq_band_dors, plot_freq_band_lat
import pickle

fnr = FileNameRetriever(bids_root)


# 1. Frequency disttribution with Minimum norm
for subj in subjects:
    if concat_raws:
        concat_file = fnr.get_filename(subj=subj, file="concat")
        filebase = concat_file.split("/")[-1].split(".")[0]
        subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
        src = fnr.get_filename(subj=subj, file=use_source_model_for_freq)
        bem_sol = fnr.get_filename(subj=subj, file="3-layer-BEM-sol")
        trans = fnr.get_trans_file(subj, concat_file)
        fwd_name = fnr.get_filename(subj=subj, file="fwd")
        freq_MNE_folder = fnr.get_filename(subj=subj, file="freqMNE")
        all_stcs_filename = (filebase + '-stc-psd-MNE.pkl')
        all_stcs_filename = os.path.join(freq_MNE_folder, all_stcs_filename)
        sensor_psd_filename = (filebase + '-sensor-psd-MNE.pkl')
        sensor_psd_filename = os.path.join(freq_MNE_folder, sensor_psd_filename)
        if not os.path.isfile(all_stcs_filename) or not os.path.isfile(sensor_psd_filename):
            raw = mne.io.read_raw(concat_file, preload=True)
            if os.path.isfile(fwd_name):
                fwd = mne.read_forward_solution(fwd_name)
            else:    
                fwd = mne.make_forward_solution(raw.info, src=src, bem=bem_sol,
                                        trans=trans, 
                                        meg=pick_meg, eeg=pick_eeg, mindist=0.2, 
                                        ignore_ref=False, 
                                        n_jobs=n_jobs, verbose=True)
                mne.write_forward_solution(fwd_name, fwd)
            noise_cov = mne.compute_raw_covariance(raw, method="empirical", n_jobs=n_jobs)
            inv = make_inverse_operator(raw.info, forward=fwd, noise_cov=noise_cov, 
                                            loose=0.2, depth=0.8)

            snr = 3.
            lambda2 = 1. / snr ** 2

            stc_psd, sensor_psd = compute_source_psd(raw, inv, lambda2=lambda2, 
                                                method='MNE', 
                                                fmin=1, fmax=45, n_fft=2048, n_jobs=n_jobs, 
                                                return_sensor=True, verbose=True)
            sensor_psd_filename = (filebase + '-sensor-psd-MNE.pkl')
            sensor_psd_filename = os.path.join(freq_MNE_folder, sensor_psd_filename)
            pickle.dump(stc_psd, open(all_stcs_filename, "wb"))
            pickle.dump(sensor_psd, open(sensor_psd_filename, "wb"))
        else:
            stc_psd = pickle.load(open(all_stcs_filename, "rb"))
            sensor_psd = pickle.load(open(sensor_psd_filename, "rb"))
        
        # Visulisation
        topos = dict()
        stcs = dict()
        topo_norm = sensor_psd.data.sum(axis=1, keepdims=True)
        stc_norm = stc_psd.sum()
        x_hemi_freq = dict
        
        for band, limits in freq_bands.items():         # normalize...
            data = sensor_psd.copy().crop(*limits).data.sum(axis=1, keepdims=True)
            topos[band] = mne.EvokedArray(100 * data / topo_norm, sensor_psd.info)
            stcs[band] = 100 * stc_psd.copy().crop(*limits).sum() / stc_norm.data

        figure = dict()
        brain = dict()
        mne.viz.set_3d_backend('mayavi')
        for band in freq_bands.keys():
            brain[band] = plot_freq_band_dors(stcs[band], band=band, subject=subj, 
                                                        subjects_dir=subjects_dir,
                                                        filebase=filebase)
            freqfilename3d = (filebase + '_' + band + '_freq_topomap_3d_dors.png')
            freqfilename3d = os.path.join(freq_MNE_folder, freqfilename3d)
            image = brain[band].save_image(freqfilename3d)
            brain[band] = plot_freq_band_lat(stcs[band], band=band, subject=subj, 
                                                        subjects_dir=subjects_dir,
                                                        filebase=filebase)
            freqfilename3d = (filebase + '_' + band + '_freq_topomap_3d_lat.png')
            freqfilename3d = os.path.join(freq_MNE_folder, freqfilename3d)
            image = brain[band].save_image(freqfilename3d)


"""
# 2. Frequency distribution with DICS beamformer:  --> to be fixed
for subj in subjects:
    if concat_raws:
        concat_file = fnr.get_filename(subj=subj, file="concat")
        filebase = concat_file.split("/")[-1].split(".")[0]
        subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
        src = fnr.get_filename(subj=subj, file=use_source_model_for_freq)
        bem_sol = fnr.get_filename(subj=subj, file="3-layer-BEM-sol")
        trans = fnr.get_trans_file(subj, concat_file)
        fwd_name = fnr.get_filename(subj=subj, file="fwd")
        freq_DICS_folder = fnr.get_filename(subj=subj, file="freqDICS")
        all_stcs_filename = (filebase + '-stc-psd-DICS.pkl')
        all_stcs_filename = os.path.join(freq_DICS_folder, all_stcs_filename)
        sensor_psd_filename = (filebase + '-sensor-psd-DICS.pkl')
        sensor_psd_filename = os.path.join(freq_DICS_folder, sensor_psd_filename)
        brain = dict()
        if not os.path.isfile(all_stcs_filename) or not os.path.isfile(sensor_psd_filename):
            raw = mne.io.read_raw(concat_file)
            print("raw loaded")
            epochs = mne.make_fixed_length_epochs(raw, duration=15, verbose=True, preload=True)
            fwd = mne.make_forward_solution(epochs.info, src=src, bem=bem_sol,
                                            trans=trans, 
                                            meg=pick_meg, eeg=pick_eeg, mindist=0.2, 
                                            ignore_ref=False, 
                                            n_jobs=n_jobs, verbose=True)
            ### convert frequency bins to logspace:
            log_freq_bands = dict()
            for band in freq_bands:
                lower, upper = freq_bands[band]
                band_range = upper - lower +1
                log_freq_bands[band] = (np.logspace(np.log10(lower), np.log10(upper), band_range))
            
            csds = dict()
            stc_DICS = dict()
            for band in freq_bands:
                csds[band] = mne.time_frequency.csd_morlet(epochs, log_freq_bands[band], 
                                                            tmin=0, tmax=15, n_jobs=n_jobs)
                csds[band] = csds[band].mean()
                # create an empty csd as csd_noise
                csd_noise = csds[band].copy()
                noise_data = np.zeros_like(csds[band].get_data())
                csd_noise.data = noise_data
                DICS = mne.beamformer.make_dics(epochs.info, fwd, csds[band],noise_csd=csd_noise, pick_ori="max-power")
                stc_DICS[band], _freq = mne.beamformer.apply_dics_csd(csds[band], DICS)
                brain[band] = plot_freq_band_dors(stc_DICS[band], band=band, subject=subj, 
                                                        subjects_dir=subjects_dir,
                                                        filebase=filebase)
                freqfilename3d = (filebase + '_' + band + '_DICS_freq_topomap_3d_dors.png')
                freqfilename3d = os.path.join(freq_DICS_folder, freqfilename3d)
                image = brain[band].save_image(freqfilename3d)
                brain[band] = plot_freq_band_lat(stc_DICS[band], band=band, subject=subj, 
                                                            subjects_dir=subjects_dir,
                                                            filebase=filebase)
                freqfilename3d = (filebase + '_' + band + '_DICS_freq_topomap_3d_lat.png')
                freqfilename3d = os.path.join(freq_DICS_folder, freqfilename3d)
                image = brain[band].save_image(freqfilename3d)
"""





"""
To do:
- make DICS Solution work :-)
- calculate diffenece between MNE and DICS solution --> is there any?
    normalize stc_psd-data for both solutions, subtract one from the other and plot
    for visual prototyping
        are there vast differences?


"""
