import pandas as pd
import os
import pickle
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.FusionNoCCDSID import FusionNoCCDSID
from pyliftover import LiftOver


def from_hg18_to_hg19(chr, coord):
    """
    object to perform hg18 --> hg19 conversion.
    ----------- REMEMBER that LIFT-OVER coordinates are 0-based!!!
    ----------- ADD +1 to obtain a values in 1-based coordinate!!
    :param chr: chromosome name, e.g. 'chr6'
    :param coord: integer, e.g. 10000
    :return: coord in hg coordinates system
    """
    lo = LiftOver('hg18', 'hg19')
    conv = lo.convert_coordinate(chr, int(coord)+1)
    hg19_coord = conv[0][1]
    return hg19_coord


def build_test_set():
    fusion_list = []
    version = 'grch37'

    #####################################
    # read ChimerDB3.0_ChimerSeq
    #####################################
    path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    file = pd.read_excel(os.path.join(path, "data/ChimerDB3.0_ChimerSeq.xlsx"),
                         index_col=0, header=0)

    # selecting only ChimerDB2.0, our test set for Onco
    file_red = file.loc[file['Source'] == 'ChimerDB 2.0']

    for i in range(len(file_red)):
        # for each fusion
        print("Test set Onco is building: %0.2f%% complete" % (i/len(file_red)*100))

        # find 5p and convert from hg18 to hg19
        chr5p = file_red.iloc[i]['H_chr']
        coord5p = from_hg18_to_hg19(chr5p, file_red.iloc[i]['H_position'])

        # find 3p and convert from hg18 to hg19
        chr3p = file_red.iloc[i]['T_chr']
        coord3p = from_hg18_to_hg19(chr3p, file_red.iloc[i]['T_position'])
        tissue = file_red.iloc[i]['Cancertype_or_disease']

        # create and save object
        fusion_list.append(FusionNoCCDSID(chr5p, coord5p, chr3p, coord3p, tissue, version))

        # TRANSLATE nucleotides sequences into amino acids sequences and add LABELS
        fusion_list[i].proteins = []
        fusion_list[i].labels = 1

        if len(fusion_list[i].sequences) != 0:
            for seq in fusion_list[i].sequences:
                sequence = seq

                # if there is a partial codon, delete it
                if len(sequence) % 3 != 0:
                    sequence = sequence[:-(len(sequence) % 3)]
                seq = Seq(sequence, generic_dna)
                prot = str(seq.translate()).split('*')[0]

                # add protein to fusion obj.proteins only if not present
                if prot not in fusion_list[i].proteins:
                    fusion_list[i].proteins.append(prot)

        with open(os.path.join(path, "data/test_set.data"), 'wb+') as f:
            pickle.dump(fusion_list, f)

    #####################################
    # add fusions in Normal Tissues, taken from TopHat-fusion applied onto Bodymap
    #####################################

    file = pd.read_excel(os.path.join(path, "data/test_set_normal_TopHat_onto_BodyMap.xlsx"),
                         header=0)

    for i in range(len(file)):
        # for each fusion
        print("Test set NotOnco from TopHatFusion is building: %0.2f%% complete" % (i/len(file)*100))

        # find 5p
        chr5p = file.iloc[i]['Chromosomes (left-right)'].split('-')[0]
        coord5p = file.iloc[i]['5p position']

        # find 3p
        chr3p = file.iloc[i]['Chromosomes (left-right)'].split('-')[1]
        coord3p = file.iloc[i]['3p position']
        tissue = file.iloc[i]['SAMPLE ID']

        # create and save object
        fusion_list.append(FusionNoCCDSID(chr5p, coord5p, chr3p, coord3p, tissue, version))

        # TRANSLATE nucleotides sequences into amino acids sequences and add LABELS
        fusion_list[i+len(file_red)].proteins = []
        fusion_list[i+len(file_red)].labels = 0

        if len(fusion_list[i+len(file_red)].sequences) != 0:
            for seq in fusion_list[i+len(file_red)].sequences:
                sequence = seq

                # if there is a partial codon, delete it
                if len(sequence) % 3 != 0:
                    sequence = sequence[:-(len(sequence) % 3)]
                seq = Seq(sequence, generic_dna)
                prot = str(seq.translate()).split('*')[0]

                # add protein to fusion obj.proteins only if not present
                if prot not in fusion_list[i+len(file_red)].proteins:
                    fusion_list[i+len(file_red)].proteins.append(prot)

        with open(os.path.join(path, "data/test_set.data"), 'wb+') as f:
            pickle.dump(fusion_list, f)

    #####################################
    # add fusions in Normal Tissues, taken from STAR-fusion applied onto Bodymap2
    #####################################
    file_1 = pd.read_csv(os.path.join(path, "data/fusions_Bodymap2.txt"), sep='\t', header=None)
    version = 'grch38'

    for i in range(len(file_1)):
        # for each fusion
        print("Test set NotOnco from STAR-Fusion is building: %0.2f%% complete" % (i/len(file_1)*100))

        # find 5p
        chr5p = file_1.iloc[i][5].split(':')[0]
        coord5p = file_1.iloc[i][5].split(':')[1]

        # find 3p
        chr3p = file_1.iloc[i][7].split(':')[0]
        coord3p = file_1.iloc[i][7].split(':')[1]
        tissue = 'None'

        # create and save object
        fusion_list.append(FusionNoCCDSID(chr5p, coord5p, chr3p, coord3p, tissue, version))

        # TRANSLATE nucleotides sequences into amino acids sequences and add LABELS
        fusion_list[i+len(file_red)+len(file)].proteins = []
        fusion_list[i+len(file_red)+len(file)].labels = 0

        if len(fusion_list[i+len(file_red)+len(file)].sequences) != 0:
            for seq in fusion_list[i+len(file_red)+len(file)].sequences:
                sequence = seq

                # if there is a partial codon, delete it
                if len(sequence) % 3 != 0:
                    sequence = sequence[:-(len(sequence) % 3)]
                seq = Seq(sequence, generic_dna)
                prot = str(seq.translate()).split('*')[0]

                # add protein to fusion obj.proteins only if not present
                if prot not in fusion_list[i+len(file_red)+len(file)].proteins:
                    fusion_list[i+len(file_red)+len(file)].proteins.append(prot)

        with open(os.path.join(path, "data/test_set.data"), 'wb+') as f:
            pickle.dump(fusion_list, f)

