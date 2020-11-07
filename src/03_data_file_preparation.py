# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from configuration import (subjects, n_jobs, bids_root, 
                            concat_raws, l_freq, h_freq, fir_design, s_freq,
                            do_filter, do_resample)
import glob
import mne
from utils.utils import FileNameRetriever , RawPreprocessor

fnr = FileNameRetriever(bids_root)
prepper = RawPreprocessor()

if concat_raws == True:
    for subj in subjects:
        raw_files = fnr.get_tsss_fifs(subj)
        print([raw_files])
        concat_name = fnr.get_filename(subj, "concat")
        if not os.path.isfile(concat_name):
            raws = []
            for idx, raw in enumerate(raw_files):
                rawname = "raw-" + str(idx)
                rawname = mne.io.read_raw(raw, preload=True)
                raws.append(rawname)

            raw_concat = mne.concatenate_raws(raws)

            if do_filter:
                raw_concat = prepper.filter_raw(raw_concat, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design, n_jobs=n_jobs)
            if do_resample:
                raw_concat = prepper.resample_raw(raw_concat, s_freq=s_freq)
            raw_concat.save(concat_name, overwrite=True)
            print(f"\nConcatenated file info:\n{raw_concat.info}")
        else:
            print("\nNot doing file concatenation, as concat-file exists...")

"""

# Combine eventfile + raw


for subj in subjects:
    raws = fnr.get_tsss_fifs(subj)
    for raw in raws:
        event_folder = fnr.get_filename(subj, "event_folder")
        try:
            event_file = prepper.get_event_file(event_folder, raw)
            event_file, event_dict = prepper.transform_eventfile(event_file)
            raw = prepper.combine_events_w_raw(raw, event_file)
            


        except Exception as e:
            print(f"\nSomething went wrong for: \n{raw}")
        



        eve = u.transform_event_file(eve_file, results_dir)
        raw = mne.io.read_raw_fif(filename).load_data()
        new_eve_filename = os.path.join(results_dir, eve_file_base + '_processed.csv')
        events = u.read_pandas_events_csv(new_eve_filename)
        event_id = dict()
        event_id = u.generate_spike_dict(events)
        print(f"\n\n\nEvents as loaded via pandas: {events}")
        events = u.drop_names_from_event_df(events)
        raw.add_events(events) #, stim_channel='STI014')
        raw = u.do_crop_and_filter(raw, results_dir, crop_tmax=None, n_jobs=n_jobs, filebase=filebase)
        raw = u.do_resample(raw, down_sfreq)
        raw, fig = u.do_ecg_reduction(raw, filename=filename, for_report_dir=for_report_dir)
        fname = (filebase + '_active_ssp_projections_ecg.png')
        fname = os.path.join(for_report_dir, fname)
        fig.savefig(fname)
        #raw = u.do_eog_reduction(raw, show_plot=show_plot, filename=filename, for_report_dir=for_report_dir)                      
        raw.save(raw_name, overwrite=True)


        if do_filter == True:
            raw = mne.io.read_raw(raw, preload=True)
            prepper.filter_raw(raw, l_freq=l_freq, h_freq=h_freq, fir_design=fir_design)
        if do_resample == True:
            prepper.resample_raw(raw, s_freq=s_freq)

"""