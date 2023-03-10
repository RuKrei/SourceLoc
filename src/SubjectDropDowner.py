import ipywidgets as widgets
import glob
import os


class SubjectDropDowner:
    """A Class to create a dropdown menu of subjects from the input folder.
    Parameters: subject_list: A list of subjects to be used in the dropdown menu.
    """

    def __init__(self, subject_list):
        self.subject_list = subject_list
        self.subject_dropdown = self.create_subject_dropdown()
        self.subject_dropdown_widget = self.create_subject_dropdown_widget()

    def create_subject_dropdown(self) -> widgets.Dropdown:
        subject_dropdown = widgets.Dropdown(
            options=self.subject_list,
            value=self.subject_list[0],
            description="Subject: ",
            disabled=False,
            layout=widgets.Layout(width="auto"),
        )
        return subject_dropdown

    def create_subject_dropdown_widget(self) -> widgets.VBox:
        subject_dropdown_widget = widgets.VBox([self.subject_dropdown])
        return subject_dropdown_widget

    def get_subject_dropdown_widget(self) -> widgets.VBox:
        return self.subject_dropdown_widget

    def get_subject_dropdown(self) -> widgets.Dropdown:
        return self.subject_dropdown


def get_subject_list(input_dir=None) -> list:
    analist = glob.glob(os.path.join(input_dir, "*"))
    analist = [os.path.basename(f) for f in analist if os.path.isdir(f)]
    fiflist = [
        os.path.basename(f)
        for f in glob.glob(os.path.join(input_dir, "*trans_tsss.fif"))
    ]
    fiflist = [f.split("_")[0] for f in fiflist]
    subjectlist = set(fiflist + analist)
    subjectlist = sorted([f for f in subjectlist])
    return subjectlist
