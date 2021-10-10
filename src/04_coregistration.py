# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
from os.path import isfile
from configuration import subjects_dir, derivatives_root, concat_raws, use_single_transfile, session
import mne
from utils.utils import FileNameRetriever
from mne_bids import BIDSPath, read_raw_bids
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sub", action="store", type=str, required=True)
args = parser.parse_args()

subj = args.sub

fnr = FileNameRetriever(derivatives_root)

subsubj = "sub-" + subj
subjects_dir = fnr.get_filename(subsubj, "subjects_dir")
bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", 
                            root=derivatives_root, processing="tsssTransEvePreproc")
print(f"\n\nThe following files with processing= \"tsssTransEvePreproc\" were found: {bids_derivatives.match()}\n\n")

if bids_derivatives.match() == []:
    bids_derivatives = BIDSPath(subject=subj, datatype="meg", session=session, task="resting", 
                            root=derivatives_root, processing="tsssTransNoEvePreproc")

#load data
# This should actually be:
#try:
#    raw = read_raw_bids(bids_path=bids_derivatives)
#except Exception as e:
#    print(e)
#    print("Looking for a proprocessed fif-file without Events...")
#    bids_derivatives.update(processing="tsssTransNoEvePreproc")

# loading manually, as BIDS query returns all kinds of things...
target_dir = os.path.join(derivatives_root, subsubj, "ses-resting", "meg", subsubj)
try:
    rawfile = glob.glob(target_dir + "*tsssTransEvePreproc_meg.fif")[0]              #### breaks if rawfile is not found!
except Exception as e:
    print(e)
# if no events
if not os.path.isfile(rawfile):
    rawfile = glob.glob(target_dir + "*tsssTransNoEvePreproc_meg.fif")[0]

if use_single_transfile == True:
    transfile = fnr.get_single_trans_file(subsubj)
    if isfile(transfile):
        print(f"Skipping coregistration, because a transfile ({transfile}) already exists")
    else:
        print(f"\n\n\n--> Transfile should be called: {transfile}\n\n\n")
        mne.gui.coregistration(subject=subsubj, subjects_dir=subjects_dir, inst=rawfile, advanced_rendering=False) # BIDS: inst=raw.filenames[0])


   

"""
This errors on WSLg with the following problem:

QStandardPaths: XDG_RUNTIME_DIR points to non-existing path '/home/idrael/.xdg', please create it with 0700 permissions.
Traceback (most recent call last):
  File "/home/idrael/git/SourceLoc/src/04_coregistration.py", line 56, in <module>
    mne.gui.coregistration(subject=subsubj, subjects_dir=subjects_dir, inst=rawfile, advanced_rendering=False) # BIDS: inst=raw.filenames[0])
  File "<decorator-gen-551>", line 24, in coregistration
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/mne/gui/__init__.py", line 179, in coregistration
    _check_backend()
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/mne/gui/_backend.py", line 35, in _check_backend
    from pyface.api import warning
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/api.py", line 13, in <module>
    from .about_dialog import AboutDialog
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/about_dialog.py", line 17, in <module>
    AboutDialog = toolkit_object("about_dialog:AboutDialog")
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/base_toolkit.py", line 152, in __call__
    module = import_module(mname, package)
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/importlib/__init__.py", line 127, in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/ui/qt4/about_dialog.py", line 23, in <module>
    from pyface.i_about_dialog import IAboutDialog, MAboutDialog
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/i_about_dialog.py", line 18, in <module>
    from pyface.image_resource import ImageResource
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/image_resource.py", line 19, in <module>
    ImageResource = toolkit_object("image_resource:ImageResource")
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/base_toolkit.py", line 152, in __call__
    module = import_module(mname, package)
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/importlib/__init__.py", line 127, in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/ui/qt4/image_resource.py", line 25, in <module>
    from pyface.i_image_resource import IImageResource, MImageResource
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/i_image_resource.py", line 14, in <module>
    from pyface.resource_manager import resource_manager
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/resource_manager.py", line 14, in <module>
    from pyface.resource.api import ResourceManager
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/resource/api.py", line 13, in <module>
    from .resource_manager import ResourceManager
  File "/home/idrael/anaconda3/envs/mne023/lib/python3.9/site-packages/pyface/resource/resource_manager.py", line 20, in <module>
    from importlib_resources import files
ModuleNotFoundError: No module named 'importlib_resources'




--> 
mkdir ~/.xdg
chmod 0700 ~/.xdg 
pip install importlib_resources

-->
leads to coreg program opening, but not displaying head shape or fiducials






# To do:
# - if "tsssTrans" in processing --> one transfile
# - if "raw" in processing --> specific transfile
#
# maybe only use processing="tsssTransEvePreproc" for frequency spectrum at first
#

"""