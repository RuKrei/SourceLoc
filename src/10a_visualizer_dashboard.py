#!/usr/bin/python
# Author: Rudi Kreidenhuber <Rudi.Kreidenhuber@gmail.com>
# License: BSD (3-clause)

import mne
import os
import glob
import pandas as pd
import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from configuration import subjects, subjects_dir, derivatives_root, n_jobs
import ipywidgets as widgets
import dash
import dash_core_components as dcc
import dash_html_components as html


os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt5'
mne.viz.set_3d_backend('pyvista')

subsubj = "sub-" + str(subjects[0])

base_dir = os.path.join(derivatives_root, subsubj)
eventfolders = glob.glob(base_dir + "/spikes/*")
eventnames = list(map(lambda e: e.split("/")[-1], eventfolders))
modalities = ["eLORETA_with_peaks", "dSPM"]
eventdict = dict()

heading = """
----
# Visualizer
----

### Use this server to: 
    * visualize .stcs
    * generate custom screenshots
    * generate custom time series

#### --> Subjectname:
"""
subject_to_ana = "##### " + str(subjects[0])

app = dash.Dash(__name__, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])
"""
app.layout = html.Div(children=[
    html.Div(className="row",
        children=[
            dcc.Markdown(heading),
            dcc.Markdown(subject_to_ana)
            ])
    ])
"""

app.layout = html.Div([
        # heading
        html.Div([
            dcc.Markdown(heading),
            dcc.Markdown(subject_to_ana),
        ], className="seven columns"),

        html.Div([
            html.H4('Top right'),
            dcc.Graph(id='g1', figure={'data': [{'y': [1, 2, 3]}]})
        ], className="three columns"),

        # body
        html.Div([
            html.H4('1. -> Select event:'),

            dcc.Graph(id='g2', figure={'data': [{'y': [1, 2, 3]}]})
        ], className="ten columns"),

        html.Div([
            html.H4('2. -> Select modality:'),
            dcc.RadioItems(
                options=[
                    {"label": "eLORETA with peaks", "value": "eLORETA"},
                    {"label": "dSPM", "value": "dSPM"},
                ]
            ),
            dcc.Graph(id='g3', figure={'data': [{'y': [1, 2, 3]}]})
        ], className="ten columns"),

        html.Div([
            html.H4('View:'),
            dcc.Graph(id='g4', figure={'data': [{'y': [1, 2, 3]}]})
        ], className="ten columns"),

        html.Div([
            html.H4('Column left'),
            dcc.Graph(id='g8', figure={'data': [{'y': [1, 2, 3]}]})
        ], className="five columns"),

        html.Div([
            html.H4('Column right'),
            dcc.Graph(id='g9', figure={'data': [{'y': [1, 2, 3]}]})
        ], className="five columns"),
    ], className="row")

for ef in eventfolders:
    e = ef.split("/")[-1]
    stcs = glob.glob(ef +"/stc_*")
    mod_and_file = dict()
    for s in stcs:
        for m in modalities:
            if m in s:
                mod_and_file[m] = s
        eventdict[e] = mod_and_file

df = pd.DataFrame.from_dict(eventdict, orient='index')

times = np.linspace(-20, 5, 6)
time_dict = {
            -20.0: 'a-m20',
            -15.0: 'b-m15',
            -10.0: 'c-m10',
            -5.0 : 'd-m5',
            -0.0 : 'e-m0',
             5.0 : 'f-p5'}

def save_ts():
    surfer_kwargs = dict(
                        hemi=hemi, subjects_dir=subjects_dir,
                        clim=dict(kind='value', lims=[low, medium, high]),
                        views=views,
                        time_unit='ms', size=(800, 800),
                        smoothing_steps=7)
    brains = dict()
    event_name = df.index[event_number]
    save_directory = os.path.join(base_dir, "spikes", event_name, "custom_time_series")
    modality = df.columns[index]

    for time in times:
        brains[time] = stc.plot(initial_time=time, **surfer_kwargs)
        brain_label = str(text_to_add + ' @ ' + str(time) + 'ms')
        brains[time].add_text(0.1, 0.9, brain_label, 'title', font_size=10)
        text_to_add = str((df.iloc[event_number][index]).split("/")[-1])
        custom_pic_name = (event_name + '-' + modality + '-at-' + time_dict[time] + '.png')
        custom_pic = os.path.join(save_directory, custom_pic_name)
        brains[time].save_image(custom_pic)
        brains[time].close()

def show_stc():
    stc = mne.read_source_estimate(df.loc[event_number.value][modality.value], subject=subsubj)

    surfer_kwargs = dict(
        hemi=hemi.value, subjects_dir=subjects_dir,
        clim=dict(kind='percent', lims=[85, 97, 99]), 
        #clim=dict(kind='value', lims=[1e-12, 4.5e-12, 8e-12]),
        views=views.value,
        time_unit='ms', size=(1100, 800),
        smoothing_steps=7)

    brain = stc.plot(initial_time=0., **surfer_kwargs)
    file_name = df.loc[event_number.value][modality.value]
    file_name = file_name.split('-rh')[0]
    file_name = file_name.split('-lh')[0]
    file_name = file_name.split("/")[-1]
    text_to_add = str(file_name)
    brain.add_text(0.1, 0.9, text_to_add, 'title', font_size=10)

if __name__ == '__main__':
    app.run_server(debug=True)