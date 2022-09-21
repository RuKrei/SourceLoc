#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from os.path import join as opj
import shutil
import glob
import platform
import argparse
import mne
from datetime import datetime
from pickle import load
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import linspace
import logging


class EpilepsyReportBuilder:
    def __init__(self, derivatives_root=None, subject=None, extras_dir=None):
        self.derivatives_root = derivatives_root
        self.subject = subject
        if not self.subject.startswith("sub-"):
            self.subject = "sub-" + self.subject
        self.extras_dir = extras_dir
        self.fbase = opj(self.derivatives_root, self.subject)
        self.fsrc = opj(self.fbase, "source_model")
        self.fanat = opj(self.fbase, "freesurfer")
        self.fprep = opj(self.fbase, "preprocessing")
        self.spikes = opj(self.fbase, "spikes")
        self.freq = opj(self.fbase, "frequency_distribution")
        self.conn = opj(self.fbase, "connectivity")
        self.ftrans = opj(self.fbase, "trans_files")
        self.freport = opj(self.fbase, "report")

    def _get_aquisition_date(self, raw):
        aquisition_date = raw.info['meas_date']
        year = aquisition_date.year
        month = aquisition_date.month
        day = aquisition_date.day
        return str(day) + '-' + str(month) + '-' + str(year)
    
    def _return_stc(self, event=None, modality="eLORETA"):
        name = os.path.basename(event)
        if modality == "eLORETA":
            return opj(event, "stc_" + self.subject.split("sub-")[-1] + "_" + name + "_eLORETA-lh.stc" )
        if modality == "dSPM":
            return opj(event, "stc_" + self.subject.split("sub-")[-1] + "_" + name + "_dSPM-stc.h5")

    def _plot_frequencies(self, band=None):
        if band == None:
            print("Expected a freq_band (str))")
            return
        freq_folder = self.freq
        matplotlib.rcParams["figure.facecolor"] = "black"
        freq_files = glob.glob(opj(freq_folder, "*_freq_topomap_3d_dors.png"))
        xhemi_files = glob.glob(opj(freq_folder, "*_x_hemi*.png"))
        for freq_file in freq_files:
            if str(band) in os.path.basename(freq_file):
                for xhemi_file in xhemi_files:
                    if str(band) in os.path.basename(xhemi_file):
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
                        return fig

    def _is_desired_event(self, de):
        is_a_fif = True if de.endswith(".fif") else False
        is_pickle = True if de.endswith(".pkl") else False
        ignorance = True if "ignore_me" in de else False
        point = True if de.startswith(".") else False
        triple_a = True if de == "AAA" else False
        return not (any([is_a_fif, is_pickle, ignorance, point, triple_a]))  # None of the above conditions should be true

    def _plot_time_course(self, event=None):
        custom_ts_folder = os.path.join(self.spikes, event, "custom_time_series")
        series = sorted(glob.glob(custom_ts_folder + "/*.png"))
        matplotlib.rcParams["figure.facecolor"] = "black"
        fig = plt.Figure()
        fig.set_figwidth(15)
        fig.set_figheight(15)
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
        fig.tight_layout()
        return fig
    
    def _plot_ECD_table(self, T1=None, drei=None, event=None):
        fig = plt.Figure(figsize=(15,25), dpi=150, facecolor="k")

        # T1 imgs
        ax1 = fig.add_subplot(5, 2, 1)
        ax1.set_title('minus 20 ms')
        ax1.set_xticks([])
        ax1.set_yticks([])
        mpimg_img = mpimg.imread(T1[4]) 
        ax1.imshow(mpimg_img)

        ax2 = fig.add_subplot(5, 2, 3)
        ax2.set_title('minus 15 ms')
        ax2.set_xticks([])
        ax2.set_yticks([])
        mpimg_img = mpimg.imread(T1[3]) 
        ax2.imshow(mpimg_img)

        ax3 = fig.add_subplot(5, 2, 5)
        ax3.set_title('minus 10 ms')
        ax3.set_xticks([])
        ax3.set_yticks([])
        mpimg_img = mpimg.imread(T1[2]) 
        ax3.imshow(mpimg_img)

        ax4 = fig.add_subplot(5, 2, 7)
        ax4.set_title('minus 5 ms')
        ax4.set_xticks([])
        ax4.set_yticks([])
        mpimg_img = mpimg.imread(T1[1]) 
        ax4.imshow(mpimg_img)

        ax5 = fig.add_subplot(5, 2, 9)
        ax5.set_title('peak')
        ax5.set_xticks([])
        ax5.set_yticks([])
        mpimg_img = mpimg.imread(T1[0]) 
        ax5.imshow(mpimg_img)

        # 3D imgs
        ax6 = fig.add_subplot(5, 2, 2)
        ax6.set_xticks([])
        ax6.set_yticks([])
        mpimg_img = mpimg.imread(drei[4]) 
        ax6.imshow(mpimg_img)

        ax7 = fig.add_subplot(5, 2, 4)
        ax7.set_xticks([])
        ax7.set_yticks([])
        mpimg_img = mpimg.imread(drei[3]) 
        ax7.imshow(mpimg_img)

        ax8 = fig.add_subplot(5, 2, 6)
        ax8.set_xticks([])
        ax8.set_yticks([])
        mpimg_img = mpimg.imread(drei[2]) 
        ax8.imshow(mpimg_img)

        ax9 = fig.add_subplot(5, 2, 8)
        ax9.set_xticks([])
        ax9.set_yticks([])
        mpimg_img = mpimg.imread(drei[1]) 
        ax9.imshow(mpimg_img)

        ax10 = fig.add_subplot(5, 2, 10)
        ax10.set_xticks([])
        ax10.set_yticks([])
        mpimg_img = mpimg.imread(drei[0]) 
        ax10.imshow(mpimg_img)

        fig.suptitle(str(event + ' - ECD'), fontsize=12)
        fig.tight_layout()
        matplotlib.rcParams["figure.facecolor"] = "black"
        return fig

    def create_report(self):
        if self.subject == None:
            raise Exception
        
        fif = opj(self.spikes, self.subject + "-epo.fif")
        if not os.path.isfile(fif):
            fif = opj(self.derivatives_root, self.subject, "ses-resting", "meg")
            fif = glob.glob(fif + "/*.fif")[0]              # won't work on windows, as \ would be needed
            fif = mne.io.read_raw_fif(fif)
        else:
            fif = mne.read_epochs(fif)
        now = str(datetime.now())
        aquisition_date = self._get_aquisition_date(fif)

        # logger
        logfile = opj(self.freport, "Reporter.log")
        logging.basicConfig(filename=logfile, filemode="w",
                            format="\n%(levelname)s --> %(message)s")
        rootlog = logging.getLogger()
        rootlog.setLevel(logging.INFO)

        rootlog.info(f"Now creating report for: {self.subject.split('sub-')[-1]} \n \
                MEG aquisition date --> {aquisition_date}\n \
                Report created --> {now}")

        try:
            title = (self.subject + ' _MEG_vom_' + aquisition_date + '_Befund-' + now)
            h5title = (self.subject + ' _MEG_vom_' + aquisition_date + '_Befund')
        except NameError as ne:
            rootlog.warning(f"Title setting with aquisition date failed: {ne}")
            title = (self.subject + ' _MEG_Befund-' + now)
            h5title = (self.subject + ' _MEG_Befund-' + now)
        
        report = mne.Report(subject=self.subject, subjects_dir=self.fanat, 
                        title=title, verbose=True)

        # Add title image
        try:
            cover_file = opj(self.extras_dir, "MEG_title.png")
            cover_title = self.subject + " MEG Befund"
            report.add_image(image=cover_file, title=cover_title)
        except FileNotFoundError as fnfe:
            rootlog.warning(f"MEG title page not found: {fnfe}")

        # Event selection --> omit events by renaming the folder/ adding a . in front
        desired_events = glob.glob(opj(self.spikes, "*"))
        rootlog.info(f"The following event-folders were found:\n{desired_events}")

        # Add topomaps
        epo_filename = opj(self.spikes, str(self.subject) + "-epo.fif")
        if os.path.isfile(epo_filename):
            concat_epochs = mne.read_epochs(epo_filename)
            noise_cov_file = opj(self.spikes, "Spikes_noise_covariance.pkl")
            times = linspace(-0.02, 0.01, 6)
            with open(noise_cov_file, 'rb') as f:
                noise_cov = load(f)
            for e in desired_events:
                event = os.path.basename(e)
                matplotlib.rcParams["figure.facecolor"] = "white"
                if self._is_desired_event(event):
                    viz_eve = concat_epochs[event].average().crop(-0.15, 0.1)         
                    fig = viz_eve.plot_joint(times=times, show=False)
                    title = str(event + " - Topomap")
                    report.add_figure(fig, title=title)

            # add stcs
                    modalities = ["eLORETA"]  # later also: "dSPM"?
                    rootlog.info(f"For event \"{event}\" the following stc-modalities were included: {modalities}.")
                    for modality in modalities:
                        try:
                            stc_file = self._return_stc(event=e, modality=modality)
                            title = str(event + " - " + modality)
                            report.add_stc(stc=stc_file, title=title,
                                        subject=self.subject, subjects_dir=self.fanat, 
                                        n_time_points=100)
                        except Exception as ex:
                            rootlog.warning(f"Couldn't include {modality} - stc to report because of: {ex}")

            # add ECD pics
                    rootlog.info(f"Now generating ECD-plot for {event}...")
                    generic_pics_folder = os.path.join(self.spikes, event, "generic_pics")
                    drei = sorted(glob.glob(generic_pics_folder + "/img_3d_ecd*.png"))
                    T1 = sorted(glob.glob(generic_pics_folder + "/img_ecd_*.png"))
                    matplotlib.rcParams["figure.facecolor"] = "black"
                    if drei != [] and T1 != []:
                        ECD_fig = self._plot_ECD_table(T1=T1, drei=drei, event=event)
                        caption = str(event + " - ECD")
                        report.add_figure(ECD_fig, title=caption, caption=caption)

            # add custom pics and custom time series
                    custom_pics_folder = os.path.join(self.spikes, event, "custom_pics")
                    custom_pics = glob.glob(custom_pics_folder + "/*.png")
                    custom_ts_folder = os.path.join(self.spikes, e, "custom_time_series")
                    custom_ts = glob.glob(custom_ts_folder + "/*.png")
                    if custom_pics is not []:
                        rootlog.info(f"Now adding custom pics for {event}...")
                        for cst in custom_pics:
                            cst_title = cst.split('/')[-1]
                            cst_title = cst_title.split('.')[0]
                            caption = event + ' - ' + cst_title
                            report.add_image(cst, title=cst_title, caption=caption)
                    if custom_ts is not []:   
                        rootlog.info(f"Now adding custom time series for {event}...") 
                        for _ in custom_ts:
                            caption = event + ' - Time course'
                            fig = plt.figure(figsize=(30, 30), dpi=150, facecolor="k")
                            fig = self._plot_time_course(event=e)
                            plt.tight_layout()
                            report.add_figure(fig, title=caption, caption=caption)
                            break

        # add frequency distribution
            #freq_file = opj(self.freq, self.subject + "_Freqs-stc-psd-MNE.pkl")   # --> would be nice, but doesn't work, freq = time-index
            #with open(freq_file, "rb") as f:
            #    stc_freqs = load(f)
            #title = str(self.subject.split("sub-")[-1] + " - Frequency distribution")
            #report.add_stc(stc=stc_freqs, title=title,  # tags=("Frequency distribution"),
            #                        subject=self.subject, subjects_dir=self.fanat)
        
        freq_bands = ["delta", "theta", "alpha", "beta", "gamma"]                #the frequency bands of interest for the analysis    
        rootlog.info(f"Now adding frequency distribution...")
        for freq in freq_bands:
            try:
                fig = self._plot_frequencies(freq)
                title = str(self.subject.split("sub-")[-1] + " - Frequency distribution - " + freq)
                report.add_figure(fig, title=title)
            except Exception as ex:
                rootlog.warning(f"Something went wrong trying to add freqs: {ex}")
        
        # BEM
        try:
            rootlog.info(f"Now adding BEM.")
            report.add_bem(self.subject, decim=4, subjects_dir=self.fanat, 
                                width=256)
        except ValueError as ve:
            rootlog.info("Could not add BEM to report, maybe a spherical model was used? Error was: {ve}")

        # Add disclaimer image
        try:
            rootlog.info(f"Now adding Disclaimer.")
            disclaimer_file = opj(self.extras_dir, 'MEG_disclaimer.png')
            report.add_image(disclaimer_file, title='disclaimer')   
        except FileNotFoundError as fnfe:
            rootlog.warning(f"Disclaimer file not found - {e}")
            
        # Save all
        try:
            rootlog.info(f"Saving...")
            title = (self.subject + " _MEG_Befund.html")
            save_name_html = os.path.join(self.freport, title)
            save_name_h5 = os.path.join(self.freport, (h5title + '.h5'))   
            report.save(save_name_html, overwrite=True)
            #report.save(save_name_h5)
            rootlog.info(f"Saving complete.")
        except Exception as ex:
            rootlog.warning(f"Saving failed because of: {ex}")

def main():
    parser = argparse.ArgumentParser()     
    parser.add_argument("--derivatives_root", action="store", type=str, required=False, 
                        help="Specify the BIDS root folder")
    parser.add_argument("--subject", action="store", type=str, required=True, 
                        help="Subject name")
    parser.add_argument("--extras_dir", action="store", type=str, required=False, 
                        help="Extras directory")

    args = parser.parse_args()

    # additional arguments
    subject = args.subject
    if not subject.startswith("sub-"):
        subject = "sub-" + subject

    if args.extras_dir:
        extras_directory = args.extras_dir
    else:
        extras_directory = os.environ.get("EXTRAS_DIRECTORY")
    
    if args.derivatives_root:
        derivatives_root = args.derivatives_root
    else:
        derivatives_root = os.environ.get("DERIVATIVES_ROOT")   
    
    reporter = EpilepsyReportBuilder(derivatives_root=derivatives_root, subject=subject, 
                                    extras_dir=extras_directory)
    reporter.create_report()
    
if __name__ == '__main__':
    main()