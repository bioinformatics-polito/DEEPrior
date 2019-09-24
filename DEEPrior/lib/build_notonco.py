import pandas as pd
import os
import pickle
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.Fusion import Fusion


def build_notonco():
    path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    file = pd.read_excel(os.path.join(path, "data/list_of_fusions_in_normal_tissues.xlsx"),
                         sheet_name='Table S2', header=1)

    fusion_list = []
    version = 'grch37'

    for i in range(len(file)):
        print("Not Onco database is building: %0.2f%% complete" % (i/len(file)*100))
        chr5p = file.loc[i]['up_chr']
        coord5p = file.loc[i]['up_Genome_pos']
        chr3p = file.loc[i]['dw_chr']
        coord3p = file.loc[i]['dw_Genome_pos']
        tissue = file.loc[i]['sample source']

        fusion_list.append(Fusion(chr5p, coord5p, chr3p, coord3p, tissue, version))

        with open(os.path.join(path, "data/NotOnco.data"), 'wb+') as f:
            pickle.dump(fusion_list, f)
