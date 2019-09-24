from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import copy
import random as rd


def sequences_by_fusion_pair(onco, not_onco_pruned):
    """
    This function groups fusion objects into dictionaries in which the fusion_pair is used as primary key of the
    dictionary.
    :param onco: list of onco objects
    :param not_onco_pruned: list of not_onco objects
    :return: onco_dict, not_onco_dict
    """
    #################################################
    # onco dict
    #################################################
    onco_dict = {}
    tot_seq = 0

    # for each onco object
    for i in onco:
        sequence = i.sequence
        if sequence != '':
            # if there is a partial codon, delete it
            if len(sequence) % 3 != 0:
                sequence = sequence[:-(len(sequence) % 3)]
            seq = Seq(sequence, generic_dna)
            prot = str(seq.translate()).split('*')[0]
            i.protein = prot

            # select only prot with len > 5!!
            if len(prot) > 5:
                # if the fusion_pair is not in the onco dictionary, add it
                if i.fusion_pair not in onco_dict.keys():
                    onco_dict[i.fusion_pair] = [prot]
                    tot_seq += 1

                # else add to the dictionary if it is not duplicated
                else:
                    if prot not in onco_dict[i.fusion_pair]:
                        onco_dict[i.fusion_pair].append(prot)
                        tot_seq += 1

    #################################################
    # not onco dict
    #################################################
    not_onco_final = copy.deepcopy(not_onco_pruned)
    not_onco_dict = {}
    el_to_pop = []
    tot_protein_sequences = 0

    ###
    # first add only protein sequences coming from recurrent fusions >= 3:
    for i in range(len(not_onco_final)):
        elem = not_onco_final[i]
        if len(elem.tissue.split('_')) > 2:
            sequences = elem.sequences
            for sequence in sequences:
                if sequence != '':
                    # if there is a partial codon, delete it
                    if len(sequence) % 3 != 0:
                        sequence = sequence[:-(len(sequence) % 3)]

                    # perform dna sequence translation into protein sequence
                    seq = Seq(sequence, generic_dna)
                    prot = str(seq.translate()).split('*')[0]
                    elem.protein = prot

                    # only if len prot is grater than 5:
                    if len(prot) > 5:
                        # if the fusion_pair is not in the not_onco dictionary, add it
                        if elem.fusion_pair not in not_onco_dict.keys():
                            not_onco_dict[elem.fusion_pair] = [prot]
                            tot_protein_sequences += 1

                        # else add to the dictionary if it is not duplicated
                        else:
                            if prot not in not_onco_dict[elem.fusion_pair]:
                                not_onco_dict[elem.fusion_pair].append(prot)
                                tot_protein_sequences += 1

            el_to_pop.append(i)

    ####
    # popping recurrent >=3 fusion objects
    for i in sorted(el_to_pop, reverse=True):
        del not_onco_final[i]

    ####
    # randomly selecting fusion pairs from not onco till the overall sample is produced
    while tot_protein_sequences < tot_seq:
        i = rd.choice(range(len(not_onco_final)))
        elem = not_onco_final[i]
        if len(elem.tissue.split('_')) > 1:
            sequences = elem.sequences
            for sequence in sequences:
                if sequence != '':
                    # if there is a partial codon, delete it
                    if len(sequence) % 3 != 0:
                        sequence = sequence[:-(len(sequence) % 3)]

                    # perform dna sequence translation into protein sequence
                    seq = Seq(sequence, generic_dna)
                    prot = str(seq.translate()).split('*')[0]
                    elem.protein = prot

                    # only if len prot is grater than 5:
                    if len(prot) > 5:

                        # if the fusion_pair is not in the not_onco dictionary, add it
                        if elem.fusion_pair not in not_onco_dict.keys():
                            not_onco_dict[elem.fusion_pair] = [prot]
                            tot_protein_sequences += 1

                        # else add to the dictionary if it is not duplicated
                        else:
                            if prot not in not_onco_dict[elem.fusion_pair]:
                                not_onco_dict[elem.fusion_pair].append(prot)
                                tot_protein_sequences += 1

        # deleting the randomly chosen object
        not_onco_final.pop(i)

    return onco_dict, not_onco_dict
