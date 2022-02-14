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
    main_protein = None
    protein_len = None
    early_stop = None
    complete_5p = None
    complete_3p = None

    a = [x.split('*')[0] for x in fusion_object.proteins_details]
    if len(a) != 0:
        ind = [len(x) for x in a].index(max([len(x) for x in a]))

        # main protein and main protein length
        main_protein = a[ind]
        protein_len = len(main_protein)

        # check if an early stop codon occurs in the protein
        if '*' in fusion_object.proteins_details[ind]:
            early_stop = 'Yes'
        else:
            early_stop = 'No'

        # check if 5p gene is complete
        if fusion_object.info_breakpoints[ind][0] == 'UTR_3p':
            complete_5p = 'Yes'
        else:
            complete_5p = 'No'

        # check if 3p gene is complete
        if fusion_object.info_breakpoints[ind][1] == 'UTR_5p':
            complete_3p = 'Yes'
        else:
            complete_3p = 'No'

    return main_protein, protein_len, early_stop, complete_5p, complete_3p


def generating_fusion_list(file, version, mode):
    fusion_list = []
    for i in range(len(file)):
        # for each fusion
        print("Dataset is building: %0.2f%% complete. Fusion number %d" % ((i + 1) / len(file) * 100, i+1))
        if str(file.iloc[i, 0]) != 'nan' and str(file.iloc[i, 1]) != 'nan' and str(file.iloc[i, 2]) != 'nan' and str(file.iloc[i, 3]) != 'nan':
            chr5p, coord5p, chr3p, coord3p = str(file.iloc[i, 0]), str(file.iloc[i, 1]), \
                                             str(file.iloc[i, 2]), str(file.iloc[i, 3])
            tissue = 'Unknown'
        # create and save object
            fusion_list.append(FusionNoCCDSID(chr5p, coord5p, chr3p, coord3p, tissue, i, version))

            # TRANSLATE nucleotides sequences into amino acids sequences and add LABELS
            fusion_list[i].proteins = []
            fusion_list[i].proteins_details = []
            fusion_list[i].main_protein = ''
            fusion_list[i].main_protein_len = None
            fusion_list[i].early_stop = None
            fusion_list[i].complete_5p = None
            fusion_list[i].complete_3p = None

            if mode == 'inference':
                fusion_list[i].labels = 'To be evaluated by the model'
            elif mode == 'retraining':
                fusion_list[i].labels = file.iloc[i, 4]

            if len(fusion_list[i].sequences) != 0:
                for seq in fusion_list[i].sequences:
                    sequence = seq

                    # if there is a partial codon, delete it
                    if len(sequence) % 3 != 0:
                        sequence = sequence[:-(len(sequence) % 3)]
                    seq = Seq(sequence, generic_dna)
                    b = str(seq.translate())
                    prot = b.split('*')[0]

                    # add protein fusion obj.proteins_details
                    fusion_list[i].proteins_details.append(b)

                    # add protein to fusion obj.proteins only if not present
                    if prot not in fusion_list[i].proteins:
                        fusion_list[i].proteins.append(prot)

                # calculating main protein, main protein length, early stopping in the protein,
                # 5p and 3p gene coming complete only if there are sequences!!
                fusion_list[i].main_protein, \
                fusion_list[i].main_protein_len, \
                fusion_list[i].early_stop, \
                fusion_list[i].complete_5p, \
                fusion_list[i].complete_3p = main_sequence(fusion_list[i])
        else:
            fusion_list.append("Missing information, fusion skipped")
            print(f'Fusion number {i+1} skipped due to lack of necessary information')

    return fusion_list


def create_fusion_list(file, tool, version, mode):
    if mode == 'inference':
        if tool == 'ChimPIPE':
            to_parse = file['juncCoord'].apply(lambda x: x.split(':'))
            five_p = to_parse.apply(lambda x: x[0])
            three_p = to_parse.apply(lambda x: x[1])
            file['chr5p'] = five_p.apply(lambda x: x.split('_')[0])
            file['coord5p'] = five_p.apply(lambda x: x.split('_')[1])
            file['chr3p'] = three_p.apply(lambda x: x.split('_')[0])
            file['coord3p'] = three_p.apply(lambda x: x.split('_')[1])
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'DeFuse':
            file['chr5p'] = 'chr' + file['gene_chromosome1'].apply(str)
            file['chr3p'] = 'chr' + file['gene_chromosome2'].apply(str)
            file['coord5p'] = file['genomic_break_pos1'].apply(int).apply(str)
            file['coord3p'] = file['genomic_break_pos2'].apply(int).apply(str)
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'EricScript':
            file['chr5p'] = 'chr' + file['chr1']
            file['chr3p'] = 'chr' + file['chr2']
            file['coord5p'] = file['Breakpoint1']
            file['coord3p'] = file['Breakpoint2']
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'FusionCatcher':
            to_parse = file['Fusion_point_for_gene_1(5end_fusion_partner)'].apply(lambda x: x.split(':'))
            file['chr5p'] = 'chr' + to_parse.apply(lambda x: x[0])
            file['coord5p'] = to_parse.apply(lambda x: x[1])

            to_parse = file['Fusion_point_for_gene_2(3end_fusion_partner)'].apply(lambda x: x.split(':'))
            file['chr3p'] = 'chr' + to_parse.apply(lambda x: x[0])
            file['coord3p'] = to_parse.apply(lambda x: x[1])
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'InFusion':
            file['chr5p'] = 'chr' + file['ref1'].apply(str)
            file['chr3p'] = 'chr' + file['ref2'].apply(str)
            file['coord5p'] = file['break_pos1']
            file['coord3p'] = file['break_pos2']
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'JAFFA':
            file['chr5p'] = file['chrom1']
            file['chr3p'] = file['chrom2']
            file['coord5p'] = file['base1']
            file['coord3p'] = file['base2']
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'SOAPfuse':
            file['chr5p'] = file['up_chr']
            file['chr3p'] = file['dw_chr']
            file['coord5p'] = file['up_Genome_pos']
            file['coord3p'] = file['dw_Genome_pos']
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'STAR-Fusion':
            file['chr5p'] = file['LeftBreakpoint'].apply(lambda x: x.split(':')[0])
            file['coord5p'] = file['LeftBreakpoint'].apply(lambda x: x.split(':')[1])

            file['chr3p'] = file['RightBreakpoint'].apply(lambda x: x.split(':')[0])
            file['coord3p'] = file['RightBreakpoint'].apply(lambda x: x.split(':')[1])

            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'TopHat':
            file['chr5p'] = file[2]
            file['chr3p'] = file[5]
            file['coord5p'] = file[3]
            file['coord3p'] = file[6]
            reduced_file = file[['chr5p', 'coord5p', 'chr3p', 'coord3p']]

        if tool == 'general':
            reduced_file = file

    if mode == 'retraining':
        reduced_file = file

    fusion_list = generating_fusion_list(reduced_file, version, mode)

    return fusion_list
