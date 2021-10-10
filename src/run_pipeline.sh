#!/bin/bash

#bst



# usage = ./run_pipeline.sh SUBJECTNAME 



python 00_prep_input.py --sub $1
python 01_create_DS_and_Folders.py --sub $1
python 02_freesurfer_and_BEM.py --sub $1
python 03_data_file_preparation.py --sub $1
python 04_coregistration.py --sub $1
python 05_frequency_spectrum.py --sub $1
#python 06_hippocampal_segmentation.py --sub $1
python 07_source_localization.py --sub $1
#python -i 08_connectivity.py --sub $1  #not implemented yet
python 09_report.py --sub $1
#
#echo "Report has been opened"
#echo "You can now create screenshots and inspect the time series via double-click on"
#echo "BIDS-Directory  --> /derivatives/SUBJECTNAME/report/10_visualizer.ipynb"
#echo
#echo
#
#echo "Finished? Update report? Then type y (y, *.&%$v.') "
#read postproc

#if [ "$postproc" = y ]; then
#    python 09_report.py
#fi



