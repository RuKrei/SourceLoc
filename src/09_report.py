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
#matplotlib.rcParams["figure.facecolor"] = "white"

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
    
    fig.faceclolor = "black"
    return fig

def plot_ECD_table(T1_imgs=None, drei_D_imgs=None, event='GR_1'):
    fig = plt.Figure(figsize=(15,25))
    
    # T1 imgs
    ax1 = fig.add_subplot(5, 2, 1)
    ax1.set_title('minus 20 ms')
    ax1.set_xticks([])
    ax1.set_yticks([])
    mpimg_img = mpimg.imread(T1_imgs[4]) 
    ax1.imshow(mpimg_img)
    
    ax2 = fig.add_subplot(5, 2, 3)
    ax2.set_title('minus 15 ms')
    ax2.set_xticks([])
    ax2.set_yticks([])
    mpimg_img = mpimg.imread(T1_imgs[3]) 
    ax2.imshow(mpimg_img)
    
    ax3 = fig.add_subplot(5, 2, 5)
    ax3.set_title('minus 10 ms')
    ax3.set_xticks([])
    ax3.set_yticks([])
    mpimg_img = mpimg.imread(T1_imgs[2]) 
    ax3.imshow(mpimg_img)
    
    ax4 = fig.add_subplot(5, 2, 7)
    ax4.set_title('minus 5 ms')
    ax4.set_xticks([])
    ax4.set_yticks([])
    mpimg_img = mpimg.imread(T1_imgs[1]) 
    ax4.imshow(mpimg_img)
    
    ax5 = fig.add_subplot(5, 2, 9)
    ax5.set_title('peak')
    ax5.set_xticks([])
    ax5.set_yticks([])
    mpimg_img = mpimg.imread(T1_imgs[0]) 
    ax5.imshow(mpimg_img)
    
    
    # 3D imgs
    ax6 = fig.add_subplot(5, 2, 2)
    #ax6.set_title('minus 20 ms')
    ax6.set_xticks([])
    ax6.set_yticks([])
    mpimg_img = mpimg.imread(drei[4]) 
    ax6.imshow(mpimg_img)
    
    ax7 = fig.add_subplot(5, 2, 4)
    #ax7.set_title('minus 15 ms')
    ax7.set_xticks([])
    ax7.set_yticks([])
    mpimg_img = mpimg.imread(drei[3]) 
    ax7.imshow(mpimg_img)
    
    ax8 = fig.add_subplot(5, 2, 6)
    #ax8.set_title('minus 10 ms')
    ax8.set_xticks([])
    ax8.set_yticks([])
    mpimg_img = mpimg.imread(drei[2]) 
    ax8.imshow(mpimg_img)
    
    ax9 = fig.add_subplot(5, 2, 8)
    #ax9.set_title('minus 5 ms')
    ax9.set_xticks([])
    ax9.set_yticks([])
    mpimg_img = mpimg.imread(drei[1]) 
    ax9.imshow(mpimg_img)
    
    ax10 = fig.add_subplot(5, 2, 10)
    #ax10.set_title('peak')
    ax10.set_xticks([])
    ax10.set_yticks([])
    mpimg_img = mpimg.imread(drei[0]) 
    ax10.imshow(mpimg_img)
    
    
    fig.suptitle(str(event + ' - ECD'), fontsize=12)
    fig.tight_layout()
    matplotlib.rcParams["figure.facecolor"] = "black"
    return fig



for subj in subjects:
    subsubj = "sub-" + subj
    meg_folder = os.path.join(derivatives_root, subsubj, "ses-resting", "meg")
    report_folder = fnr.get_filename(subsubj, "report")
    spike_folder = fnr.get_filename(subsubj, "spikes")
    fif = glob.glob(meg_folder + "/*finalEpochs_meg.fif")
    print(f" fif --> {fif}")
    raw = mne.io.read_raw(fif[-1])
    aquisition_date = raw.info['meas_date']
    year = aquisition_date.year
    month = aquisition_date.month
    day = aquisition_date.day
    aquisition_date = (str(day) + '-' + str(month) + '-' + str(year))

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
    #bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", 
    #                    root=derivatives_root, processing="finalEpochs_meg")    
    #raw = read_raw_bids(bids_derivatives)
    events, event_ids = mne.events_from_annotations(raw)
    print (f"event_ids: {event_ids}")
    to_drop = []
    for e in event_ids:
        if e.lower() == "ignore_me" or e.upper() == "AAA" or e.startswith("."):
            to_drop.append(e)
    if len(to_drop) > 0:
        for todr in to_drop:
            del event_ids[todr]
    print (f"event_ids (after deletion of unwanted events): {event_ids}")
    epochs = mne.Epochs(raw, events=events, event_id=event_ids, tmin=-1.5, tmax=1, baseline=(-1.5,-1), 
                        on_missing = "ignore", event_repeated="merge")
    events = epochs.event_id
    
    if events.keys() != []:
        spike_folder = fnr.get_filename(subsubj, "spikes")
        desired_events = glob.glob(spike_folder + "/*")
        #print (desired_events)
        """    
        for e in desired_events:
            e = e.split("/")[-1]
        """
        for e in events:
            if e.lower() == "ignore_me" or e.upper() == "AAA" or e.startswith("."): # pointless double check
                print (f"Omitting {e} from Analysis")
            elif e in events:
                print(f"Adding data from {e} to report...")
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
                # add ECD-picks
                drei = sorted(glob.glob(generic_pics_folder + "/img_3d_ecd*.png"))
                T1 = sorted(glob.glob(generic_pics_folder + "/img_ecd_*.png"))
                ECD_fig = plot_ECD_table(T1_imgs=T1, drei_D_imgs=drei, event='e')
                
                # custom pics
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
                # ECD fig
                caption = str(e) + ' --> Equivalent current dipole model'
                report.add_figs_to_section(ECD_fig, captions=caption, section=e)
                # custom pics
                if custom_pics is not []:
                    for cst in custom_pics:
                        cst_title = cst.split('/')[-1]
                        cst_title = cst_title.split('.')[0]
                        caption = e + ' --> ' + cst_title
                        report.add_images_to_section(cst, section=e, captions=caption)
                #matplotlib.rcParams["figure.facecolor"] = "black"
                if custom_ts is not []:    
                    for cts in custom_ts:
                        caption = e + ' --> Time course'
                        fig = plt.figure(figsize=(30, 30), dpi=150, facecolor="k")
                        fig = plot_time_course(sorted(custom_ts), event=e)
                        plt.tight_layout()
                        report.add_figs_to_section(fig, section=e, captions=caption)
                        break
                                                                                                                                 
    
    freq_folder = fnr.get_filename(subsubj, "freqMNE")
    matplotlib.rcParams["figure.facecolor"] = "black"
    for band in freq_bands.keys():
        freq_files = glob.glob(freq_folder + '/*_freq_topomap_3d_dors.png')
        xhemi_files = glob.glob(freq_folder + '/*_x_hemi*.png')
        for freq_file in freq_files:
            if str(band) in freq_file.split("/")[-1]:
                for xhemi_file in xhemi_files:
                    if str(band) in xhemi_file.split("/")[-1]:
                        fig = plt.Figure(facecolor="k")
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

                        matplotlib.rcParams["figure.facecolor"] = "black"
                        fig.tight_layout()
                        

                        # add to report
                        cap = 'Freq.-Distribution --> ' + str(band)
                        sec = "Freqs"
                        report.add_figs_to_section(fig, section=sec, captions=cap)


    # Addendum
    sec="Addendum"
    # Preprocessing data
    preproc_folder = fnr.get_filename(subsubj, "preprocessing")
    addendum_files = glob.glob(preproc_folder + "/*.png")
    
    for f in sorted(addendum_files):
        if f != []:
            cap = f.split("/")[-1]
            report.add_images_to_section(f, captions=cap, section=sec)
    
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
postprocessor in order to append custom pics etc.
explain in report why something is analysed + clinical relevance
    i.e. unilateral focal slowing in MEG 95-100% chance of ipsilateral seizure onset...
    add relevant literature
add disclaimer + title file to report dir (in 01_create_DS_and_folders)
"""