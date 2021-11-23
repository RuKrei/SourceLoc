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

    def create_report(self):
        if self.subject == None:
            raise Exception
        
        fif = opj(self.spikes, self.subject + "-epo.fif")
        fif = mne.read_epochs(fif)
        now = str(datetime.now())
        aquisition_date = self._get_aquisition_date(fif)

        try:
            title = (self.subject + ' _MEG_vom_' + aquisition_date + '_Befund-' + now)
            h5title = (self.subject + ' _MEG_vom_' + aquisition_date + '_Befund')
        except NameError:
            print("Title setting with date failed")
            title = (self.subject + ' _MEG_Befund-' + now)
            h5title = (self.subject + ' _MEG_Befund-' + now)
        
        print(f"\n\n\n\n\n\nNow creating report for: {self.subject.split('sub-')[-1]} \n \
                MEG aquisition date --> {aquisition_date}\n \
                Report created --> {now}")
        
        report = mne.Report(subject=self.subject, subjects_dir=self.fanat, 
                        title=title, verbose=True)


    # Add title image
        try:
            cover_file = opj(self.extras_dir, "MEG_title.png")
            cover_title = self.subject + " MEG Befund"
            report.add_images_to_section(cover_file, section=cover_title, captions=cover_title)
        except FileNotFoundError as e:
            print(e)

    # Event selection --> omit events by renaming the folder/ adding a . in front
        desired_events = glob.glob(opj(self.spikes, "*"))
        print(f"Desired events are: {desired_events}")

    # Add topomaps 

    # To do:    add topomaps from cropped stcs, save (cropped) evokeds in spike-subfolder
    #           use those as stcs for eLORETA below



    
        epo_filename = opj(self.spikes, str(self.subject) + "-epo.fif")
        concat_epochs = mne.read_epochs(epo_filename)
        noise_cov_file = opj(self.spikes, "Spikes_noise_covariance.pkl")
        times = linspace(-0.02, 0.01, 6)
        with open(noise_cov_file, 'rb') as f:
            noise_cov = load(f)
        for de in desired_events:
            de = os.path.basename(de)
            if self._is_desired_event(de):
                viz_eve = concat_epochs[de].average().crop(-0.15, 0.1)         
                fig = viz_eve.plot_joint(times=times, show=False)
                title = str(de + " - Topomap")
                report.add_figure(fig, title=title)

    # add stcs
        for e in desired_events:
            modalities = ["eLORETA"]  # later also: "dSPM"?
            event = os.path.basename(e)
            for modality in modalities:
                try:
                    stc_file = self._return_stc(event=e, modality=modality)
                    title = str(self.subject.split("sub-")[-1] + " - " + modality + " - " + event)
                    report.add_stc(stc=stc_file, title=title,
                                subject=self.subject, subjects_dir=self.fanat, 
                                n_time_points=100)
                except Exception as e:
                    print(e)
    
    # add frequency distribution
        #freq_file = opj(self.freq, self.subject + "_Freqs-stc-psd-MNE.pkl")   # --> would be nice, but doesn't work, freq = time-index
        #with open(freq_file, "rb") as f:
        #    stc_freqs = load(f)
        #title = str(self.subject.split("sub-")[-1] + " - Frequency distribution")
        #report.add_stc(stc=stc_freqs, title=title,  # tags=("Frequency distribution"),
        #                        subject=self.subject, subjects_dir=self.fanat)
        
        freq_bands = ["delta", "theta", "alpha", "beta", "gamma"]                #the frequency bands of interest for the analysis    
        for freq in freq_bands:
            try:
                fig = self._plot_frequencies(freq)
                title = str(self.subject.split("sub-")[-1] + " - Frequency distribution - " + freq)
                report.add_figure(fig, title=title)
            except Exception as e:
                print(e)
        
    # BEM
        try:
            report.add_bem_to_section(self.subject, decim=4, subjects_dir=self.fanat, 
                                    section='BEM')
        except ValueError:
            print ("Could not add BEM to report, it seems a spherical model was used...")

    # Add disclaimer image
        try:
            disclaimer_file = opj(self.extras_dir, 'MEG_disclaimer.png')
            report.add_images_to_section(disclaimer_file, section='disclaimer', captions='End notes')   
        except FileNotFoundError as e:
            print(e)
            
    # Save all
        title = (self.subject + " _MEG_Befund.html")
        save_name_html = os.path.join(self.freport, title)
        save_name_h5 = os.path.join(self.freport, (h5title + '.h5'))   
        report.save(save_name_html, overwrite=True)
        #report.save(save_name_h5)








def main():
    parser = argparse.ArgumentParser()     
    parser.add_argument("--derivatives_root", action="store", type=str, required=False, 
                        help="Specify the BIDS root folder")
    parser.add_argument("--subject", action="store", type=str, required=True, 
                        help="Subject name")
    parser.add_argument("--extras_dir", action="store", type=str, required=True, 
                        help="Extras directory")

    args = parser.parse_args()
    
    reporter = EpilepsyReportBuilder(derivatives_root=args.derivatives_root, subject=args.subject, 
                                    extras_dir=args.extras_dir)
    reporter.create_report()

    print("\n\n\nall good")
    

if __name__ == '__main__':
    main()