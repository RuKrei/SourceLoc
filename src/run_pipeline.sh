#!/bin/bash

#bst

#python /home/idrael/DATA/git/SourceLoc/src/00_prep_input.py
#python /home/idrael/DATA/git/SourceLoc/src/01_create_DS_and_Folders.py
#python /home/idrael/DATA/git/SourceLoc/src/02_freesurfer_and_BEM.py
#python /home/idrael/DATA/git/SourceLoc/src/03_data_file_preparation.py
#python /home/idrael/DATA/git/SourceLoc/src/04_coregistration.py
#python /home/idrael/DATA/git/SourceLoc/src/05_frequency_spectrum.py
#python /home/idrael/DATA/git/SourceLoc/src/06_hippocampal_segmentation.py
#python /home/idrael/DATA/git/SourceLoc/src/07_source_localization.py
# python -i /home/idrael/DATA/git/SourceLoc/src/08_connectivity.py
#python /home/idrael/DATA/git/SourceLoc/src/09_report.py

echo "You can now create screenshots and inspect the time series via double-click on"
echo "BIDS-Directory  --> /derivatives/SUBJECTNAME/report/10_visualizer.ipynb"
echo
echo

echo "Finished? Update report? Then type y (y, *.&%$v.') "
read postproc

if [ "$postproc" = y ]; then
    python /home/idrael/DATA/git/SourceLoc/src/09_report.py
fi