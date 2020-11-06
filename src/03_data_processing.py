# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from configuration import subjects, openmp, n_jobs, bids_root, concat_raws
import glob
import mne
from utils.utils import FileNameRetriever

fnr = FileNameRetriever(bids_root)

if concat_raws == True:
    for subj in subjects:
        raws = fnr.get_tsss_fifs(subj)

        rawfiles = []
        for idx, raw in enumerate(raws):
            rawname = "raw" + str(idx)
            rawname = mne.io.read_raw(raw)
            rawfiles.append(rawname)

        raw_concat = mne.concatenate_raws(rawfiles)
        print(raw_concat.info)