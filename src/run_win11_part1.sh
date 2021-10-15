#!/bin/bash

# usage = ./run_pipeline.sh SUBJECTNAME 

python 00_prep_input.py --sub $1
python 01_create_DS_and_Folders.py --sub $1
python 02_freesurfer_and_BEM.py --sub $1
python 03_data_file_preparation.py --sub $1