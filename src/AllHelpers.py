
def convert_dcm_folder(subj):
    try:
        anafolder = os.path.join(data_root, subj, subj)
        folder = str(glob.glob((anafolder + '/1*/100*/100*'), recursive=True)[0])
        convert_directory(folder, anafolder, compression=True, reorient=True)
    except Exception as e:
        print(e)