# SourceLoc

This is still under construction. Once finished it should become a software pipeline to create and maintain a BIDS dataset of 
clinical MEG resting state measurements and for electromagnetic source imaging in python with free software. 

Under the hood it mainly uses 
 - mne-python
 - mne-bids
 - nipype
 - nilearn

## Installation
 - clone the git repository (git clone https://github.com/RuKrei/SourceLoc.git)
 - install mne (https://mne.tools/dev/install/index.html)
 - Activate mne-environment (conda activate name_of_your_mne_environment)
 - pip install -r requirements.txt

## Configuration
 - freesurfer (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) should be installed on your system for use of anatomical processing steps. The pipeline checks, if SUBJECTS_DIR is set.

 - BIDSROOT (root directory for output of the pipeline) and INPUTFOLDER (where MRI and MEEG-data are stored) need to be set. This can be done either when launching a process as arguments, or in form of environment variables.

 - OPENMP (the number of processes to use) defaults to 1 and should be set to meet your needs.

 - SRCSPACING (https://mne.tools/dev/overview/cookbook.html#setting-up-source-space) should be set on your system, it defaults to "oct6".

## Running the pipeline
 - put MEEG-data in INPUT_FOLDER, filenames MUST contain the name of the subject
 - put MRI-data in a folder named with subjects name in INPUT_FOLDER
 - cd to directory of pipeline on harddisk
 - python run_pipeline.py
 - - if no arguments are given you will be prompted for subject name to process 
 - - specify INPUT_FOLDER and BIDS_ROOT via one of the options below

## BIDS_ROOT and INPUT_FOLDER (mandatory) +
## OPENMP and SRCSPACING (optional)
### permanent solution --> set shell variables (shown for bash or zsh)
 - export BIDS_ROOT=/path/to/store/results/in
 - export INPUT_FOLDER=/path/to/input/folder
### temporary solution --> provide as arguments
 - python run_pipeline.py --bidsroot /path/to/store/results/in --inputfolder /path/to/input/folder --openmp 8 --srcspacing ico4