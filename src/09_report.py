#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

from mne import Report
import mne
import os
import glob
import numpy as np
from datetime import datetime
import pdb
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from configuration import (subjects, n_jobs, bids_root, freq_bands)
from utils.utils import FileNameRetriever

fnr = FileNameRetriever(bids_root)

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
    subjects_dir = fnr.get_filename(subj=subj, file="subjects_dir")
    meg_folder = fnr.get_filename(subj, "meg")
    report_folder = fnr.get_filename(subj, "report")
    epos = glob.glob(meg_folder + "/*epo.fif")
    fifs = glob.glob(meg_folder + "/*task-resting*.fif")
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

    report = Report(subject=subj, subjects_dir=subjects_dir, 
                        title=title, verbose=True, raw_psd=True)
    
    # Add title image
    cover_file = report_folder + '/MEG_title.png'
    cover_title = (subj + ' MEG Befund')
    report.add_images_to_section(cover_file, section=cover_title, captions=cover_title) 

    for epo in epos:
        epochs = mne.read_epochs(epo)
        events = epochs.event_id
        if events.keys() != []:
            for e in events.keys():
                if e == "ignore_me" or e == "AAA" or e.startswith("."):
                    pass
                else:
                    viz_eve = epochs[e].average().crop(-0.2, 0.2)
                    times = np.linspace(-0.02, 0.01, 6)
                    figs = viz_eve.plot_joint(times=times, show=False)
                    for f in figs:
                        cap = e + ' --> Topomaps'
                        report.add_figs_to_section(f, captions=cap, section=e)
    

    # Frequenzverteilung
    freq_folder = fnr.get_filename(subj, "freqMNE")
    for band in freq_bands.keys():
        freq_files = glob.glob(freq_folder + '/*_freq_topomap_3d_dors.png')
        xhemi_files = glob.glob(freq_folder + '/*_x_hemi*.png')
        for freq_file in freq_files:
            if str(band) in freq_file.split("/")[-1]:
                for xhemi_file in xhemi_files:
                    if str(band) in xhemi_file.split("/")[-1]:
                        fig = plt.Figure()
                        #fig.suptitle((' Frequenzverteilung - ' + str(band)), fontsize=12)
                        fig.set_figwidth(15)
                        fig.set_figheight(15)
                        
                        # dorsal
                        ax1 = fig.add_subplot(2, 2, 3)
                        ax1.set_title('dorsal')
                        ax1.set_xticks([])
                        ax1.set_yticks([])
                        mpimg_img = mpimg.imread(freq_file) 
                        ax1.imshow(mpimg_img)

                        # lateral
                        ax2 = fig.add_subplot(2, 2, 1)
                        ax2.set_title('left hemisphere')
                        ax2.set_xticks([])
                        ax2.set_yticks([])
                        lat_file = (freq_file.split("dors.png")[0] + "lat_lh.png")
                        mpimg_img = mpimg.imread(lat_file) 
                        ax2.imshow(mpimg_img)

                        ax3 = fig.add_subplot(2, 2, 2)
                        ax3.set_title('right hemisphere')
                        ax3.set_xticks([])
                        ax3.set_yticks([])
                        lat_file = (freq_file.split("dors.png")[0] + "lat_rh.png")
                        mpimg_img = mpimg.imread(lat_file) 
                        ax3.imshow(mpimg_img)

                        # xhemi
                        ax4 = fig.add_subplot(2, 2, 4)
                        ax4.set_title('xhemi')
                        ax4.set_xticks([])
                        ax4.set_yticks([])
                        mpimg_img = mpimg.imread(xhemi_file) 
                        ax4.imshow(mpimg_img)

                        fig.tight_layout()

                        # add to report
                        cap = 'Freq.-Distribution --> ' + str(band)
                        sec = "Freqs"
                        report.add_figs_to_section(fig, section=sec, captions=cap)



                        ### read images, arrang in figure, give title, add to report
                

    # Add disclaimer image
    disclaimer_file = report_folder + '/MEG_disclaimer.png'
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