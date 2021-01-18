# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import numpy as np
import matplotlib.pyplot as plt
from configuration import (subjects, n_jobs, derivatives_root, use_source_model_for_sourceloc, inv_loose_option,
                            use_fwd_model_for_sourceloc, do_volume_source_loc, source_loc_methods,
                            pick_meg, pick_eeg, concat_raws, signal_to_noise_ratio, minimum_norm_ori,
                            peaks_tmin, peaks_tmax, peaks_mode, peaks_nr_of_points, dip_times, session,
                            derivatives_root, subjects_dir)
import mne
from mne_bids import BIDSPath, read_raw_bids
from utils.utils import FileNameRetriever, RawPreprocessor, get_peak_points
import glob
from nilearn.plotting import plot_anat

mne.viz.set_3d_backend("pyvista")

prepper = RawPreprocessor()
fnr =FileNameRetriever(derivatives_root)

snr = signal_to_noise_ratio
lambda2 = 1. / snr ** 2

for subj in subjects:
    subsubj = "sub-" + subj
    bem_sol = fnr.get_filename(subj=subsubj, file=use_fwd_model_for_sourceloc)
    trans_file = fnr.get_single_trans_file(subsubj)

    bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", root=derivatives_root, processing="tsssTransEvePreproc")
    print(f"\nUsing the following file for source localization: {bids_derivatives.match()}")
    raw = read_raw_bids(bids_derivatives)
    events, event_ids = mne.events_from_annotations(raw)
    epochs = mne.Epochs(raw, events=events, event_id=event_ids, tmin=-1.5, tmax=1, baseline=(-1.5,-1), on_missing = "ignore")

    noise_cov = mne.compute_covariance(epochs, tmax=-1, 
                                    #projs=, 
                                    method='auto',
                                    n_jobs=n_jobs)

    data_cov = mne.compute_covariance(epochs,
                                    tmin=-0.5, 
                                    tmax=0.3, 
                                    #projs=, 
                                    method='auto',
                                    n_jobs=n_jobs)
    
    for event in event_ids.keys():
        eventname = str(event)
        if eventname == "ignore_me" or eventname == "AAA":
            pass
        else:
            e = epochs[eventname].load_data().average()
            spike_folder = fnr.get_filename(subsubj, "spikes")
            e_folder = os.path.join(spike_folder, eventname)
            cp_folder = os.path.join(spike_folder, eventname, "custom_pics")
            cts_folder = os.path.join(spike_folder, eventname, "custom_time_series")
            gp_folder = os.path.join(spike_folder, eventname, "generic_pics")
            folders = [e_folder, cp_folder, cts_folder, gp_folder]
            if not os.path.isdir(e_folder):
                for f in folders:
                    print(f)
                    os.mkdir(f)
            
            src_file = fnr.get_filename(subsubj, use_source_model_for_sourceloc)
            if os.path.isfile(src_file):
                src = mne.read_source_spaces(src_file)
            else:
                print("Source model not found, aborting...")
            print(bem_sol)
            fwd = mne.make_forward_solution(e.info, src=src, bem=bem_sol,
                                        trans=trans_file, 
                                        meg=pick_meg, eeg=pick_eeg, mindist=0.2, 
                                        ignore_ref=False, 
                                        n_jobs=n_jobs, verbose=True)
                         
            inv = mne.minimum_norm.make_inverse_operator(e.info, forward=fwd, noise_cov=noise_cov, 
                                        loose=inv_loose_option, depth=0.8)
    

            # Distributed source models

            for m in source_loc_methods:
                stc_name = "stc_" + m + "_" + eventname
                stc_name = 'stc_' + m


                if m == 'dSPM':
                    stc_name = mne.minimum_norm.apply_inverse(e, inv, lambda2,
                                method='dSPM', pick_ori='vector')
                    surfer_kwargs = dict(hemi='split', subjects_dir=subjects_dir,
                                clim=dict(kind='percent', lims=[90, 96, 99.85]),
                                views=['lat', 'med'], 
                                colorbar=True,
                                initial_time=0, time_unit='ms', 
                                size=(1000, 800), smoothing_steps=10)
                    
                    brain = stc_name.plot(**surfer_kwargs)
                    label = str(subj + " - " + eventname + " - Vector solution")
                    brain.add_text(0.1, 0.9, label, 'title', font_size=10)
                    
                    img_f_name = ('img_stc_' + subj + '_' + eventname + '_' + m + '.png')
                    img_f_name = os.path.join(gp_folder, img_f_name)
                    brain.save_image(img_f_name)
                    stc_f_name = ('stc_' + subj + '_' + eventname + '_' + m)
                    stc_f_name = os.path.join(e_folder, stc_f_name)
                    e.save(stc_f_name)
                

                else:
                    stc_name = mne.minimum_norm.apply_inverse(e, inv, lambda2,
                                method=m, pick_ori=minimum_norm_ori)
                    surfer_kwargs = dict(hemi='split', subjects_dir=subjects_dir,
                                clim=dict(kind='percent', lims=[90, 96, 99.85]),
                                views=['lat', 'med'], 
                                colorbar=True,
                                initial_time=0, time_unit='ms', 
                                size=(1000, 800), smoothing_steps=10)
                    
                    brain = stc_name.plot(**surfer_kwargs)
                    label = str(subj + " - " + eventname + " - " +  m)
                    brain.add_text(0.1, 0.9, label, 'title', font_size=10)
                    
                    img_f_name = ('img_stc_' + subj + '_' + eventname + '_' + m + '.png')
                    img_f_name = os.path.join(gp_folder, img_f_name)
                    brain.save_image(img_f_name)
                    stc_f_name = ('stc_' + subj + '_' + eventname + '_' + m)
                    stc_f_name = os.path.join(e_folder, stc_f_name)
                    e.save(stc_f_name)
                    
                    if m == "eLORETA":
                        rh_peaks = get_peak_points(stc_name, hemi='rh', tmin=peaks_tmin, 
                                                    tmax=peaks_tmax, nr_points=peaks_nr_of_points, mode=peaks_mode)
                        lh_peaks = get_peak_points(stc_name, hemi='lh', tmin=peaks_tmin, 
                                                    tmax=peaks_tmax, nr_points=peaks_nr_of_points, mode=peaks_mode)
                        label = str(subj + " - " + eventname + " - " +  m + " - max. activation points")
                        brain.add_text(0.1, 0.9, label, 'title', font_size=10)
                        for p in rh_peaks:
                            brain.add_foci(p, color='green', coords_as_verts=True, hemi='rh', scale_factor=0.6, alpha=0.9)
                        for p in lh_peaks:
                            brain.add_foci(p, color='green', coords_as_verts=True, hemi='lh', scale_factor=0.6, alpha=0.9)
                        stc_f_name = ('stc_' + subj + '_' + eventname + '_' + m + "_with_peaks")
                        stc_f_name = os.path.join(e_folder, stc_f_name)
                        stc_name.save(stc_f_name)
                        img_f_name = ('img_stc_' + subj + '_' + eventname + '_' + m + '_with_peaks.png')
                        img_f_name = os.path.join(gp_folder, img_f_name)
                        brain.save_image(img_f_name)

            # Dipoles
            for start, stop in dip_times.values():
                dip_epoch = e.copy().crop(start, stop).pick('meg')
                ecd = mne.fit_dipole(dip_epoch, noise_cov, bem_sol, trans=trans_file)[0]
                best_idx = np.argmax(ecd.gof)
                best_time = ecd.times[best_idx]
                trans = mne.read_trans(trans_file)
                mri_pos = mne.head_to_mri(ecd.pos, mri_head_t=trans, subject=subsubj, subjects_dir=subjects_dir)
                t1_file_name = os.path.join(subjects_dir, subsubj, 'mri', 'T1.mgz')
                stoptime = str(abs(int(stop*1000)))
                if stoptime == "5":
                    stoptime = "05"
                title = str(eventname + ' - ECD @ minus ' + stoptime + ' ms')
                t1_fig = plot_anat(t1_file_name, cut_coords=mri_pos[0], title=title)
                t1_f_name_pic = ('img_ecd_' + eventname + '_' + '_Dipol_' + stoptime + '.png')
                t1_f_name_pic = os.path.join(e_folder, "generic_pics", t1_f_name_pic)
                t1_fig.savefig(t1_f_name_pic)
                fig_3d = ecd.plot_locations(trans, subsubj, subjects_dir, mode="orthoview")
                fig_3d_pic = ('img_3d_ecd_' + eventname + '_' + '_Dipol_' + stoptime + '.png')
                fig_3d_pic = os.path.join(e_folder, "generic_pics", fig_3d_pic)
                fig_3d.savefig(fig_3d_pic)
                plt.close("all")




"""
To do:

- volume source localization
- Plot difference between prediction and result: https://mne.tools/stable/auto_tutorials/source-modeling/plot_dipole_fit.html?highlight=ecd
"""
