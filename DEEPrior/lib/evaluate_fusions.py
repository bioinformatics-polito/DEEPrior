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
                      '3pCommonName', '3pEnsg', '3pGeneFunctionality', '3pGeneDescription', 'MainProtein', 'Proteins']

    # 1. load the tokenizer
    tokenizer = load_tokenizer()

    # 2. load the model
    print("Loading the model..")
    model = load_model(trained_model_path)

    df_res = pd.DataFrame(columns=df_res_columns)
    y_test = []
    i = 0

    # 3. for each fusion in the list, extract the protein sequences
    for fusion in list_fusions:
        i = i+1
        print("Deep learning processing: %0.2f%% completed" %(i/len(list_fusions)*100))
        x_test, sequences, label = retrieve_test_sequences(fusion, tokenizer)
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
                                    'MainProteins': main_protein_seq,
                                    'Proteins': sequences
                                    }, ignore_index=True)
        else:
            try:
                df_res = df_res.append({'FusionPair': fusion.fusion_pair,
                                        'OncogenicProbability': 'Not Applicable, not protein coding genes',
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
                                        'MainProteins': main_protein_seq,
                                        'Proteins': sequences
                                        }, ignore_index=True)
            except IndexError:
                pass
            continue
    # 6. save the results to file
    df_res.to_csv(output_file_path, sep='\t')
    print("Results saved at: "+output_file_path)


