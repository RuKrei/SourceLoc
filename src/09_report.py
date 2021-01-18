#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from mne import Report
import mne
from mne_bids import read_raw_bids, BIDSPath
import os
import glob
import numpy as np
from datetime import datetime
import pdb
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from configuration import (subjects, subjects_dir, n_jobs, derivatives_root, extras_dir, freq_bands, session)
from utils.utils import FileNameRetriever

mne.viz.set_3d_backend("pyvista")
matplotlib.rcParams["figure.facecolor"] = "black"

fnr = FileNameRetriever(derivatives_root)

def plot_time_course(series, event='GR_1', filename=None):
    fig.suptitle(str(event + ' Time course'), fontsize=12)
    
    ax1 = fig.add_subplot(3, 2, 1)
    ax1.set_title('minus 20 ms')
    ax1.set_xticks([])
    ax1.set_yticks([])
    mpimg_img = mpimg.imread(sorted(series)[0]) 
    ax1.imshow(mpimg_img)
    
    ax2 = fig.add_subplot(3, 2, 2)
    ax2.set_title('minus 15 ms')
    ax2.set_xticks([])
    ax2.set_yticks([])
    mpimg_img = mpimg.imread(sorted(series)[1]) 
    ax2.imshow(mpimg_img)
    
    ax3 = fig.add_subplot(3, 2, 3)
    ax3.set_title('minus 10 ms')
    ax3.set_xticks([])
    ax3.set_yticks([])
    mpimg_img = mpimg.imread(sorted(series)[2]) 
    ax3.imshow(mpimg_img)
    
    ax4 = fig.add_subplot(3, 2, 4)
    ax4.set_title('minus 5 ms')
    ax4.set_xticks([])
    ax4.set_yticks([])
    mpimg_img = mpimg.imread(sorted(series)[3]) 
    ax4.imshow(mpimg_img)
    
    ax5 = fig.add_subplot(3, 2, 5)
    ax5.set_title('peak')
    ax5.set_xticks([])
    ax5.set_yticks([])
    mpimg_img = mpimg.imread(sorted(series)[4]) 
    ax5.imshow(mpimg_img)
    
    ax6 = fig.add_subplot(3, 2, 6)
    ax6.set_title('plus 5 ms')
    ax6.set_xticks([])
    ax6.set_yticks([])
    mpimg_img = mpimg.imread(sorted(series)[5]) 
    ax6.imshow(mpimg_img)
    
    return fig



for subj in subjects:
    subsubj = "sub-" + subj
    meg_folder = os.path.join(derivatives_root, subsubj, "ses-resting", "meg")
    report_folder = fnr.get_filename(subsubj, "report")
    spike_folder = fnr.get_filename(subsubj, "spikes")
    fifs = glob.glob(meg_folder + "/*EvePreproc*.fif")
    for fif in fifs:
        raw = mne.io.read_raw_fif(fif)
        aquisition_date = raw.info['meas_date']
        year = aquisition_date.year
        month = aquisition_date.month
        day = aquisition_date.day
        aquisition_date = (str(day) + '-' + str(month) + '-' + str(year))
        del(raw)
        break
    now = str(datetime.now())

    try:
        title = (subj + ' _MEG_vom_' + aquisition_date + '_Befund-' + now)
        h5title = title = (subj + ' _MEG_vom_' + aquisition_date + '_Befund')
    except NameError as the_error:
        print("Title setting with date failed")
        title = (subj + ' _MEG_Befund-' + now)
        h5title = (subj + ' _MEG_Befund-' + now)

    report = Report(subject=subsubj, subjects_dir=subjects_dir, 
                        title=title, verbose=True, raw_psd=True)


    # Add title image
    cover_file = extras_dir + '/MEG_title.png'
    cover_title = (subj + ' MEG Befund')
    report.add_images_to_section(cover_file, section=cover_title, captions=cover_title) 


    # Epochs
    bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", root=derivatives_root, processing="tsssTransEvePreproc")    
    raw = read_raw_bids(bids_derivatives)
    events, event_ids = mne.events_from_annotations(raw)
    epochs = mne.Epochs(raw, events=events, event_id=event_ids, tmin=-1.5, tmax=1, baseline=(-1.5,-1), on_missing = "ignore")
    events = epochs.event_id
    print("\n Events are: ", events)
    
    if events.keys() != []:
        for e in events.keys():
            spike_folder = fnr.get_filename(subsubj, "spikes")
            if e == "ignore_me" or e == "AAA" or e.startswith("."):
                pass
            else:
                # Visualize Topomaps
                cap = str(e) + " --> Topomaps"
                viz_eve = epochs[e].average().crop(-0.2, 0.2)
                times = np.linspace(-0.02, 0.01, 6)
                figs = viz_eve.plot_joint(times=times, show=False)
                
                # Add to report
                for fig in figs:
                    report.add_figs_to_section(fig, captions=cap, section=e)

                # Find generic and custom pics
                generic_pics_folder = os.path.join(spike_folder, e, "generic_pics")
                dSPM_file = glob.glob(generic_pics_folder + "/*_dSPM.png")
                #eLO_file = glob.glob(generic_pics_folder + "/*_eLORETA.png")
                eLO_peak_file = glob.glob(generic_pics_folder + "/*_eLORETA_with_peaks.png")

                custom_pics_folder = os.path.join(spike_folder, e, "custom_pics")
                custom_pics = glob.glob(custom_pics_folder + "/*.png")
                custom_ts_folder = os.path.join(spike_folder, e, "custom_time_series")
                custom_ts = glob.glob(custom_ts_folder + "/*.png")

                # Add to report
                # generics
                caption = e + ' --> dSPM'
                report.add_images_to_section(dSPM_file, captions=caption, section=e)
                caption = str(e) + ' --> eLORETA + peaks'
                report.add_images_to_section(eLO_peak_file, captions=caption, section=e)
                # custom pics
                if custom_pics is not []:
                    for cst in custom_pics:
                        cst_title = cst.split('/')[-1]
                        cst_title = cst_title.split('.')[0]
                        caption = e + ' --> ' + cst_title
                        report.add_images_to_section(cst, section=e, captions=caption)
                if custom_ts is not []:    
                    for cts in custom_ts:
                        caption = e + ' --> Time course'
                        fig = plt.figure(figsize=(30, 30), dpi=150)
                        fig = plot_time_course(sorted(custom_ts), event=e)
                        plt.tight_layout()
                        report.add_figs_to_section(fig, section=e, captions=caption)
                        break


    # Frequenzverteilung
    freq_folder = fnr.get_filename(subsubj, "freqMNE")
    for band in freq_bands.keys():
        freq_files = glob.glob(freq_folder + '/*_freq_topomap_3d_dors.png')
        xhemi_files = glob.glob(freq_folder + '/*_x_hemi*.png')
        for freq_file in freq_files:
            if str(band) in freq_file.split("/")[-1]:
                for xhemi_file in xhemi_files:
                    if str(band) in xhemi_file.split("/")[-1]:
                        fig = plt.Figure(facecolor="black")
                        #fig.suptitle((' Frequenzverteilung - ' + str(band)), fontsize=12)
                        fig.set_figwidth(15)
                        fig.set_figheight(15)
                        
                        # lateral
                        ax2 = fig.add_subplot(3, 2, 1)
                        ax2.set_title('Left hemisphere', color=(1,1,1))
                        ax2.set_xticks([])
                        ax2.set_yticks([])
                        lat_file = (freq_file.split("dors.png")[0] + "lat_lh.png")
                        mpimg_img = mpimg.imread(lat_file) 
                        ax2.imshow(mpimg_img)

                        ax3 = fig.add_subplot(3, 2, 2)
                        ax3.set_title('Right hemisphere', color=(1,1,1))
                        ax3.set_xticks([])
                        ax3.set_yticks([])
                        lat_file = (freq_file.split("dors.png")[0] + "lat_rh.png")
                        mpimg_img = mpimg.imread(lat_file) 
                        ax3.imshow(mpimg_img)

                        # medial
                        ax5 = fig.add_subplot(3, 2, 3)
                        ax5.set_title('Left hemisphere', color=(1,1,1))
                        ax5.set_xticks([])
                        ax5.set_yticks([])
                        med_file = (freq_file.split("dors.png")[0] + "med_lh.png")
                        print(med_file)
                        mpimg_img = mpimg.imread(med_file)
                        ax5.imshow(mpimg_img)

                        ax6 = fig.add_subplot(3, 2, 4)
                        ax6.set_title('Right hemisphere', color=(1,1,1))
                        ax6.set_xticks([])
                        ax6.set_yticks([])
                        med_file = (freq_file.split("dors.png")[0] + "med_rh.png")
                        mpimg_img = mpimg.imread(med_file) 
                        ax6.imshow(mpimg_img)

                         # dorsal
                        ax1 = fig.add_subplot(3, 2, 5)
                        ax1.set_title('Dorsal view', color=(1,1,1))
                        ax1.set_xticks([])
                        ax1.set_yticks([])
                        mpimg_img = mpimg.imread(freq_file) 
                        ax1.imshow(mpimg_img)

                        # xhemi
                        ax4 = fig.add_subplot(3, 2, 6)
                        ax4.set_title('Cross hemisphere comparison', color=(1,1,1))
                        ax4.set_xticks([])
                        ax4.set_yticks([])
                        mpimg_img = mpimg.imread(xhemi_file) 
                        ax4.imshow(mpimg_img)

                        fig.tight_layout()

                        # add to report
                        cap = 'Freq.-Distribution --> ' + str(band)
                        sec = "Freqs"
                        report.add_figs_to_section(fig, section=sec, captions=cap)


    # Preprocessing data
    preproc_folder = fnr.get_filename(subsubj, "preprocessing")
    EOG_files = glob.glob(preproc_folder + "/*EOG*")
    ECG_files = glob.glob(preproc_folder + "/*ECG*")
    coreg_files = glob.glob(preproc_folder + "/*coreg*")
    
    try:
        report.add_bem_to_section(subsubj, decim=4, n_jobs=n_jobs, subjects_dir=subjects_dir, section='BEM')
    except ValueError:
        print ("Could not add BEM to report, it seems a spherical model was used...")
                

    # Add disclaimer image
    disclaimer_file = extras_dir + '/MEG_disclaimer.png'
    report.add_images_to_section(disclaimer_file, section='disclaimer', captions='End notes')   

    ### Save all
    save_name_html = os.path.join(report_folder, (title + '.html'))
    save_name_h5 = os.path.join(report_folder, (h5title + '.h5'))   
    report.save(save_name_html)
    #report.save(save_name_h5)



"""
To do:
create report in .html and in .pdf format
postprocessor in order to append custom pics etc.
make beautiful + informative visualizations/graphs
explain in report why something is analysed + clinical relevance
    i.e. unilateral focal slowing in MEG 95-100% chance of ipsilateral seizure onset...
    add relevant literature

add disclaimer + title file to report dir (in 01_create_DS_and_folders)
"""