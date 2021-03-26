#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import os
import subprocess
from configuration import subjects, derivatives_root
from utils.utils import FileNameRetriever


fnr = FileNameRetriever(derivatives_root)

def run_shell_command(command):
    subprocess.run(command, shell=True, capture_output=True, check=True)


for subj in subjects:
    subsubj = "sub-" + subj
    report_folder = fnr.get_filename(subsubj, "report")
    reporter = os.path.join(report_folder, "report.ipynb")
    
    command = str("jupyter nbconvert --to python " + reporter)
    run_shell_command(command)

    runme = os.path.join(report_folder, "report.py")
    
    command = str("python " + runme)
    run_shell_command(command)