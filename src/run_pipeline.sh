#!/bin/bash

#bst

python 01_preparations.py
python 02_freesurfer_and_BEM.py
#python 03_data_file_preparation.py
#python 04_coregistration.py
#python 05_frequency_spectrum.py
#python 06_hippocampal_segmentation.py  # not implemented yet
#python 07_source_localization.py
#python -i 08_connectivity.py   #not implemented yet

#python 09_report.py
#
#echo "Report has been opened"
#echo "You can now create screenshots and inspect the time series via double-click on"
#echo "BIDS-Directory  --> /derivatives/SUBJECTNAME/report/10_visualizer.ipynb"
#echo
#echo
#
echo "Finished? Update report? Then type y (y, *.&%$v.') "
read postproc

if [ "$postproc" = y ]; then
    python 09_report.py
fi