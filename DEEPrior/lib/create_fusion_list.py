import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.FusionNoCCDSID import FusionNoCCDSID
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def main_sequence(fusion_object):
    """
    This function takes only the most complete transcript for each fusion
    :param fusion_object: each single fusion object
    :return: main protein
    """
    protein = []
    if len(fusion_object.portions) == 2:
        if len(fusion_object.portions[0].sequences) != 0:
            if len(fusion_object.portions[1].sequences) != 0:
                list1 = fusion_object.portions[0].sequences
                list2 = fusion_object.portions[1].sequences
                sequence = max(list1, key=len) + max(list2, key=len)

                # translate into protein
                if len(sequence) % 3 != 0:
                    sequence = sequence[:-(len(sequence) % 3)]
                seq = Seq(sequence, generic_dna)
                protein.append(str(seq.translate()).split('*')[0])

    return protein


def generating_fusion_list(file, version, mode):
    fusion_list = []

    for i in range(len(file)):
        # for each fusion
        print("Dataset is building: %0.2f%% complete" % ((i+1)/len(file)*100))

        # find 5p and 3p
        chr5p, coord5p, chr3p, coord3p = str(file.iloc[i,0]), str(file.iloc[i,1]), str(file.iloc[i,2]), str(file.iloc[i,3])

        tissue = 'Unknown'

        # create and save object
        fusion_list.append(FusionNoCCDSID(chr5p, coord5p, chr3p, coord3p, tissue, version))

        # TRANSLATE nucleotides sequences into amino acids sequences and add LABELS
        fusion_list[i].main_protein = main_sequence(fusion_list[i])
        fusion_list[i].proteins = []
        if mode == 'inference':
            fusion_list[i].labels = 'To be evalated by the model'
        elif mode == 'retraining':
            fusion_list[i].labels = file.iloc[i,4]

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

    return fusion_list


def create_fusion_list(file, tool, version, mode):
    
    if mode == 'inference':
        if tool == 'DeFuse':
            file['chr5p'] = 'chr' + file['gene_chromosome1']
            file['chr3p'] = 'chr' + file['gene_chromosome2']
            file['coord5p'] = file['genomic_break_pos1'].apply(int).apply(str)
            file['coord3p'] = file['genomic_break_pos2'].apply(int).apply(str)
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]
        	
        if tool == 'STAR-Fusion':
            file['chr5p'] = file['LeftBreakpoint'].apply(lambda x: x.split(':')[0])
            file['coord5p'] = file['LeftBreakpoint'].apply(lambda x: x.split(':')[1])

            file['chr3p'] = file['RightBreakpoint'].apply(lambda x: x.split(':')[0])
            file['coord3p'] = file['RightBreakpoint'].apply(lambda x: x.split(':')[1])

            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'general':
            reduced_file = file

    if mode == 'retraining':
        reduced_file = file


    fusion_list = generating_fusion_list(reduced_file, version, mode)

    return fusion_list
