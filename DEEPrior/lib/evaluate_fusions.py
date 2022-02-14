import os
from pathlib import Path
import numpy as np
import pandas as pd
from keras.models import load_model
from lib.deep_utils import load_tokenizer,retrieve_test_sequences

MIN_PROT_SEQ = 5
MAX_PROT_SEQ = 4000
dir_path = Path(__file__).absolute().parent.parent.as_posix()


def evaluate_fusions(list_fusions, trained_model_path, output_file_path):
    df_res_columns = ['FusionPair', 'OncogenicProbability', 'Version', 'Chr5p', 'Coord5p', '5pStrand', '5pCommonName',
                      '5pEnsg', '5pGeneFunctionality', '5pGeneDescription', 'Chr3p', 'Coord3p', '3pStrand',
                      '3pCommonName', '3pEnsg', '3pGeneFunctionality', '3pGeneDescription', 'MainProteinLength',
                      'TruncatedProtein', '5p_gene_complete', '3p_gene_complete', 'MainProtein']

    # 1. load the tokenizer
    tokenizer = load_tokenizer()

    # 2. load the model
    print("Loading the model..")
    model = load_model(trained_model_path)

    df_res = pd.DataFrame(columns=df_res_columns)
    y_test = []
    i = 0
    # 3. for each fusion in the list, extract the protein sequences
    toskip = [f for f, fusion in enumerate(list_fusions) if (fusion == 'Missing information, fusion skipped' or (fusion.chr5p == 'Invalid value' or fusion.chr3p == 'Invalid value'))]
    for fusion in list_fusions:
        i = i+1
        print("Deep learning processing: %0.2f%% completed" % (i/len(list_fusions)*100))
        if (i-1) in toskip:
            df_res.loc[df_res.shape[0]] = [None] * len(df_res_columns)
            if fusion == 'Missing information, fusion skipped':
                df_res.FusionPair[i-1] = 'Missing values'
                df_res.OncogenicProbability[i-1] = 'fusion skipped'
            else:
                df_res.FusionPair[i - 1] = 'Invalid input data'
                df_res.OncogenicProbability[i - 1] = 'fusion skipped'
        else:
            x_test, sequences, label = retrieve_test_sequences(fusion, tokenizer)
            if label == -1:
                df_res = df_res.append({'FusionPair': 'No protein produced',
                                        'OncogenicProbability': 'Not Applicable',
                                        'Version': fusion.version,
                                        'Chr5p': fusion.chr5p,
                                        'Coord5p': fusion.coord5p,
                                        '5pStrand': None,
                                        '5pCommonName': None,
                                        '5pEnsg': None,
                                        '5pGeneFunctionality': None,
                                        '5pGeneDescription': None,
                                        'Chr3p': fusion.chr3p,
                                        'Coord3p': fusion.coord3p,
                                        '3pStrand': None,
                                        '3pCommonName': None,
                                        '3pEnsg': None,
                                        '3pGeneFunctionality': None,
                                        '3pGeneDescription': None,
                                        'MainProteinLength': None,
                                        'TruncatedProtein': None,
                                        '5p_gene_complete': None,
                                        '3p_gene_complete': None,
                                        'MainProtein': None
                                        }, ignore_index=True)
                continue
            else:
                try:
                    if fusion.portions[0].genes[0].strand == 1:
                        fiveprime_strand = '+'
                    elif fusion.portions[0].genes[0].strand == -1:
                        fiveprime_strand = '-'
                except IndexError:
                    fiveprime_strand = ''
                    pass
                try:
                    if fusion.portions[1].genes[0].strand == 1:
                        threeprime_strand = '+'
                    elif fusion.portions[1].genes[0].strand == -1:
                        threeprime_strand = '-'
                except IndexError:
                    threeprime_strand = ''

                try:
                    main_protein_seq = fusion.main_protein
                except (IndexError, AttributeError):
                    main_protein_seq = ''

            # 4. predict the oncogenic probability for all the protein sequences in a fusion and extract the max
            if len(x_test) > 0:
                y_test.append(label)
                labels_pred = model.predict(x_test)
                predictions = [pred[0].astype(np.float) for pred in labels_pred]

                # 5. save the results on the final file
                df_res = df_res.append({'FusionPair': fusion.fusion_pair,
                                        'OncogenicProbability': max(predictions),
                                        'Version': fusion.version,
                                        'Chr5p': fusion.chr5p,
                                        'Coord5p': fusion.coord5p,
                                        '5pStrand': fiveprime_strand,
                                        '5pCommonName': fusion.portions[0].genes[0].common_name,
                                        '5pEnsg': fusion.portions[0].ensg,
                                        '5pGeneFunctionality': fusion.protein_cod[0],
                                        '5pGeneDescription': fusion.portions[0].genes[0].description,
                                        'Chr3p': fusion.chr3p,
                                        'Coord3p': fusion.coord3p,
                                        '3pStrand': threeprime_strand,
                                        '3pCommonName': fusion.portions[1].genes[0].common_name,
                                        '3pEnsg': fusion.portions[1].ensg,
                                        '3pGeneFunctionality': fusion.protein_cod[1],
                                        '3pGeneDescription': fusion.portions[1].genes[0].description,
                                        'MainProteinLength': fusion.main_protein_len,
                                        'TruncatedProtein': fusion.early_stop,
                                        '5p_gene_complete': fusion.complete_5p,
                                        '3p_gene_complete': fusion.complete_3p,
                                        'MainProtein': main_protein_seq
                                        # 'Proteins': sequences
                                        }, ignore_index=True)
            # else:
            #     try:
            #         df_res = df_res.append({'FusionPair': fusion.fusion_pair,
            #                                 'OncogenicProbability': 'Not Applicable',
            #                                 'Version': fusion.version,
            #                                 'Chr5p': fusion.chr5p,
            #                                 'Coord5p': fusion.coord5p,
            #                                 '5pStrand': fiveprime_strand,
            #                                 '5pCommonName': fusion.portions[0].genes[0].common_name,
            #                                 '5pEnsg': fusion.portions[0].ensg,
            #                                 '5pGeneFunctionality': fusion.protein_cod[0],
            #                                 '5pGeneDescription': fusion.portions[0].genes[0].description,
            #                                 'Chr3p': fusion.chr3p,
            #                                 'Coord3p': fusion.coord3p,
            #                                 '3pStrand': threeprime_strand,
            #                                 '3pCommonName': fusion.portions[1].genes[0].common_name,
            #                                 '3pEnsg': fusion.portions[1].ensg,
            #                                 '3pGeneFunctionality': fusion.protein_cod[1],
            #                                 '3pGeneDescription': fusion.portions[1].genes[0].description,
            #                                 'MainProteinLength': fusion.main_protein_len,
            #                                 'TruncatedProtein': fusion.early_stop,
            #                                 '5p_gene_complete': fusion.complete_5p,
            #                                 '3p_gene_complete': fusion.complete_3p,
            #                                 'MainProtein': main_protein_seq
            #                                 }, ignore_index=True)
            #     except IndexError:
            #         df_res = df_res.append({'FusionPair': None,
            #                                 'OncogenicProbability': 'Not Applicable',
            #                                 'Version': fusion.version,
            #                                 'Chr5p': fusion.chr5p,
            #                                 'Coord5p': fusion.coord5p,
            #                                 '5pStrand': None,
            #                                 '5pCommonName': None,
            #                                 '5pEnsg': None,
            #                                 '5pGeneFunctionality': None,
            #                                 '5pGeneDescription': None,
            #                                 'Chr3p': fusion.chr3p,
            #                                 'Coord3p': fusion.coord3p,
            #                                 '3pStrand': None,
            #                                 '3pCommonName': None,
            #                                 '3pEnsg': None,
            #                                 '3pGeneFunctionality': None,
            #                                 '3pGeneDescription': None,
            #                                 'MainProteinLength': None,
            #                                 'TruncatedProtein': None,
            #                                 '5p_gene_complete': None,
            #                                 '3p_gene_complete': None,
            #                                 'MainProtein': None
            #                                 }, ignore_index=True)
            #
            #     continue
    # 6. save the results to file
    df_res.to_csv(output_file_path, sep='\t')
    print("Results saved at: "+output_file_path)


