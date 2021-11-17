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
        cover_file = opj(self.extras_dir, "MEG_title.png")
        cover_title = self.subject + " MEG Befund"
        report.add_images_to_section(cover_file, section=cover_title, captions=cover_title)

    # BEM
        try:
            report.add_bem_to_section(self.subject, decim=4, subjects_dir=self.fanat, section='BEM')
        except ValueError:
            print ("Could not add BEM to report, it seems a spherical model was used...")
                

    # Add disclaimer image
        disclaimer_file = opj(self.extras_dir, 'MEG_disclaimer.png')
        report.add_images_to_section(disclaimer_file, section='disclaimer', captions='End notes')   

    # Save all
        save_name_html = os.path.join(self.freport, (title + '.html'))
        save_name_h5 = os.path.join(self.freport, (h5title + '.h5'))   
        report.save(save_name_html)
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