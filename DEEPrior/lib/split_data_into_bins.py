import numpy as np
import os
import pandas as pd
import pickle


def split_data_into_bins(onco_dict, not_onco_dict, n):
    """
        This function divides Onco and NotOnco sequences in N equal bins, assuring
        equal number of Onco and NotOnco sequences in each bin and absolutely
        independence of each bin from the others on Fusion-pair name

    """
    path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))

    fusion_pair_onco = []
    fusion_pair_not_onco = []
    lens_onco = []
    lens_not_onco = []

    for i in onco_dict.keys():
        fusion_pair_onco.append(i)
        lens_onco.append(len(onco_dict[i]))

    for i in not_onco_dict.keys():
        fusion_pair_not_onco.append(i)
        lens_not_onco.append(len(not_onco_dict[i]))

    max_onco_for_bin = np.ceil(sum(lens_onco) / n)
    max_not_onco_for_bin = np.ceil(sum(lens_not_onco) / n)

    # for each bin
    for i in range(1, n + 1):
        # initialize the bin
        globals()['bin_%s' % i] = pd.DataFrame(columns=['Sequences', 'Label'])
        actual_onco = 0
        actual_not_onco = 0

        # until the number of onco sequences is lower than the expected number,
        # add onco sequences to the bin
        while actual_onco < max_onco_for_bin and len(lens_onco) > 0:
            # select the fusion_pair with the maximum number of sequences
            copy = lens_onco.copy()
            copy.sort(reverse=True)
            for n_seq in copy:
                # n_seq = max(lens_onco)
                if actual_onco + n_seq <= max_onco_for_bin:
                    fusion_pair = fusion_pair_onco[lens_onco.index(n_seq)]
                    seqs = onco_dict[fusion_pair]
                    # obj_ = FusionPair(fusion_pair, seqs, 1)
                    # globals()['bin_%s' % i].append(obj_)
                    row = pd.DataFrame([[seqs, 1]], columns=['Sequences', 'Label'], index=[fusion_pair])
                    globals()['bin_%s' % i] = globals()['bin_%s' % i].append(row)
                    actual_onco += n_seq
                    # deleting the fusion pair from the list
                    fusion_pair_onco.pop(lens_onco.index(n_seq))
                    lens_onco.pop(lens_onco.index(n_seq))

        # until the number of onco sequences is lower than the expected number,
        # add onco sequences to the bin
        while actual_not_onco < max_not_onco_for_bin and len(lens_not_onco) > 0:
            # select the fusion_pair with the maximum number of sequences
            copy = lens_not_onco.copy()
            copy.sort(reverse=True)
            for n_seq in copy:
                # n_seq = max(lens_onco)
                if actual_not_onco + n_seq <= max_not_onco_for_bin:
                    fusion_pair = fusion_pair_not_onco[lens_not_onco.index(n_seq)]
                    seqs = not_onco_dict[fusion_pair]
                    # obj_ = FusionPair(fusion_pair, seqs, 0)
                    # globals()['bin_%s' % i].append(obj_)
                    row = pd.DataFrame([[seqs, 0]], columns=['Sequences', 'Label'], index=[fusion_pair])
                    globals()['bin_%s' % i] = globals()['bin_%s' % i].append(row)
                    actual_not_onco += n_seq
                    # deleting the fusion pair from the list
                    fusion_pair_not_onco.pop(lens_not_onco.index(n_seq))
                    lens_not_onco.pop(lens_not_onco.index(n_seq))

        with open(os.path.join(path, 'data/BINS/bin_%s.pickle' % i), 'wb') as f:
            pickle.dump(globals()['bin_%s' % i], f)

        o = 0
        n = 0
        for k in range(len(globals()['bin_%s' % i])):
            if globals()['bin_%s' % i]['Label'][k] == 0:
                n += len(globals()['bin_%s' % i]['Sequences'][k])
            if globals()['bin_%s' % i]['Label'][k] == 1:
                o += len(globals()['bin_%s' % i]['Sequences'][k])

        print("Bin_%s   Onco: %s, NotOnco: %s" % (i, o, n))

    return
