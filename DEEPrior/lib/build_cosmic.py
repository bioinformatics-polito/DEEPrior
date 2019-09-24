"""
This function is intended to parse CosmicFusionExport.tsv file and to output a list of Cosmic objects written in
data/Cosmic.data
"""

import pandas as pd
import pickle
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.Cosmic import Cosmic


def build_cosmic():
    # absolute path of deepfusion_analysis folder
    path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))

    file = pd.read_csv(os.path.join(path, "data/CosmicFusionExport.tsv"), sep='\t', header=0)

    file1 = file.drop_duplicates(subset='Translocation Name').iloc[1:, :]

    # create a list of COSMIC objects for each valid translocation name and save it
    transl_list = []

    for i in range(len(file1)):
        if '(' not in file1.iloc[i]['Translocation Name'] and '?' not in file1.iloc[i]['Translocation Name']:
            print("Cosmic database is building: %0.2f%% complete" % (i/len(file1)*100))
            transl_list.append(Cosmic(file1.iloc[i]['Translocation Name'], file1.iloc[i]['Primary site'], 'grch37'))

        with open(os.path.join(path, "data/Cosmic.data"), 'wb+') as f1:
            pickle.dump(transl_list, f1)

    return
