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
# Report and visualizer should run from subjects report folder either as .ipynb (generic) or 
# bendable to the users will --> the 09_run_reoprt.py and 10_visualizer.py files are merely 
# fancy link generators

python /home/idrael/DATA/git/SourceLoc/src/09_run_report.py
#code /home/idrael/DATA/git/SourceLoc/src/10_visualizer.ipynb
#jupyter nbconvert --to python /home/idrael/DATA/git/SourceLoc/src/report.ipynb  # should be run from results dir 
#python /home/idrael/DATA/git/SourceLoc/src/report.py
#   ---> 09_report.py should execute this as a shell command
