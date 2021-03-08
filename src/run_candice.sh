#!/bin/bash

#bst
#python 01_create_DS_and_Folders.py
#python 02_freesurfer_and_BEM.py
python 03_data_file_preparation.py
python 04_coregistration.py
python 05_frequency_spectrum.py
#python 06_hippocampal_segmentation.py
python 07_source_localization.py
#python -i 08_connectivity.py
python 09_report.py
#code 10_visualizer.ipynb
