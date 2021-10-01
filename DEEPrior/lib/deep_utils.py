import os
from pathlib import Path
import configparser
import pickle
import pandas as pd
import numpy as np
import re
import tensorflow as tf
from keras.models import Sequential
from keras.preprocessing.sequence import pad_sequences
from keras.layers import MaxPooling1D, Dense, CuDNNLSTM, Dropout, Embedding, Bidirectional, Conv1D, LSTM

ONCO = 1
NOT_ONCO = 0
MIN_PROT_SEQ = 5
MAX_PROT_SEQ = 4000
NUM_AA = 20
dir_path = Path(__file__).absolute().parent.parent.as_posix()
config_path = os.path.join(dir_path, "resources/config.txt")
training_set_path = os.path.join(dir_path, "data/training_set.csv")
tokenizer_path = os.path.join(dir_path, "resources/tokenizer.pickle")


def load_config(config_section):
    """This function loads a section of the configuration file (specified in the config_section variable)"""
    print(config_section)
    config = configparser.RawConfigParser()
    config.read(config_path)
    retraining_params = dict(config.items(config_section))
    return retraining_params


def load_tokenizer():
    """This function loads the tokenizer saved at tokenizer_path."""
    with open(tokenizer_path, 'rb') as handle:
        tokenizer = pickle.load(handle)
    return tokenizer


def get_session_gpu(per_process_gpu_memory_fraction, allow_growth, visible_device_list):
    """This function  sets the GPU configuration."""
    gpu_options = tf.compat.v1.GPUOptions(per_process_gpu_memory_fraction=per_process_gpu_memory_fraction,
                                          allow_growth=allow_growth,
                                          visible_device_list=visible_device_list)
    return tf.compat.v1.Session(config=tf.compat.v1.ConfigProto(gpu_options=gpu_options))


def get_session_cpu(intra_op_parallelism_threads, inter_op_parallelism_threads):
    """This function  sets the CPU configuration."""
    session_conf = tf.compat.v1.ConfigProto(
        intra_op_parallelism_threads=intra_op_parallelism_threads,
        inter_op_parallelism_threads=inter_op_parallelism_threads
    )
    return tf.compat.v1.Session(config=session_conf)


def retrieve_test_sequences(fusion, tokenizer):
    """This function extracts the protein sequences from the fusion object and pads them.
    The function returns the list of those sequences and their labels."""
    x_test = []
    sequences = []
    label = -1
    try:
        if len(fusion.proteins) > 0:
            sequences = []
            for seq in fusion.proteins:
                if len(seq) > MIN_PROT_SEQ and len(seq) <= MAX_PROT_SEQ:
                    sequences.append(seq)
                if len(sequences) > 0:
                    x_test_nopad = tokenizer.texts_to_sequences(sequences)
                    x_test = pad_sequences(x_test_nopad, maxlen=MAX_PROT_SEQ, padding='post')
                    label = fusion.labels
    except:
        pass
    return x_test, sequences, label


def create_val_train_sets(val_percentage, fusions_list, tokenizer, add_default):
    """This function creates a validation and a training sets starting from a list of fusions.
    If the add_default flag is set to True, the new datasets are built also by using the
    default DEEPrior training set."""
    df_columns = ['FusionPair', 'ProteinSequence', 'Label', 'SequenceLength']
    train_df = pd.DataFrame(columns=df_columns)
    val_df = pd.DataFrame(columns=df_columns)

    # 1. add the fusions in the default dataset if needed
    if add_default is True:
        default_training_set = pd.read_csv(training_set_path, sep="\t")
        for i in range(len(default_training_set)):
            protein_list = re.sub("[\'\[\]]", "", default_training_set['Proteins'][i]).split(",")
            for protein_seq in protein_list:
                train_df = train_df.append({'FusionPair': default_training_set['FusionPair'][i],
                                            'ProteinSequence': protein_seq,
                                            'Label': default_training_set['Label'][i],
                                            'SequenceLength': len(protein_seq)},
                                           ignore_index=True)

    # 2. add to the new dataset the fusions retrieved from the input file
    if fusions_list is not None:
        for fusion in fusions_list:
            for protein_seq in fusion.proteins:
                if len(protein_seq) > MIN_PROT_SEQ:
                    train_df = train_df.append({'FusionPair': fusion.fusion_pair,
                                                'ProteinSequence': protein_seq,
                                                'Label': fusion.labels,
                                                'SequenceLength': len(protein_seq)},
                                               ignore_index=True)
    # 3. shuffle the dataframe
    train_df = train_df.sample(frac=1)

    # 4. compute how many sequences you need for the validation set
    len_val_set = round(len(train_df) * float(val_percentage))

    # 5. divide the data in training and validation
    val_df = train_df[:len_val_set]
    train_df = train_df[len_val_set:]

    # 6. pad the training and the validation data
    x_val_nopad = tokenizer.texts_to_sequences(np.array(val_df['ProteinSequence'].values))
    x_val = pad_sequences(x_val_nopad, maxlen=MAX_PROT_SEQ, padding='post')
    y_val = np.asarray(val_df['Label'].values).astype('int32')
    x_train_nopad = tokenizer.texts_to_sequences(np.array(train_df['ProteinSequence'].values))
    x_train = pad_sequences(x_train_nopad, maxlen=MAX_PROT_SEQ, padding='post')
    y_train = np.asarray(train_df['Label'].values).astype('int32')

    print("--- Number of validation sequences: {0} (onco: {1} and not onco: {2} )".format(str(len(val_df)), str(
        len(val_df[val_df['Label'] == ONCO])), str(len(val_df[val_df['Label'] == NOT_ONCO]))))
    print("--- Number of training sequences: {0} (onco: {1} and not onco: {2} )".format(str(len(train_df)), str(
        len(train_df[train_df['Label'] == ONCO])), str(len(train_df[train_df['Label'] == NOT_ONCO]))))

    return x_train, y_train, x_val, y_val


def cnn_bilstm_dropout2_cpu(embedding_size, embeddings_initializer, filters, kernel_size,
                            activation_function_conv1d, dropout, max_pool_win, hidden_nodes,
                            activation_function_output, loss, optimizer):
    """This function creates and compile a CNN-Bidirectional(LSTM) model with 2 levels of dropout.
    This model s on CPU only."""
    print('Creating a CNN-Bidirectional(LSTM) model with two layers of dropout (CPU)...')
    model = Sequential()
    model.add(Embedding(NUM_AA + 1, embedding_size, input_length=MAX_PROT_SEQ,
                        embeddings_initializer=embeddings_initializer))
    model.add(Conv1D(filters,
                     kernel_size,
                     padding='valid',
                     activation=activation_function_conv1d,
                     strides=1))
    model.add(MaxPooling1D(pool_size=max_pool_win))
    model.add(Dropout(dropout))
    model.add(Bidirectional(LSTM(hidden_nodes)))
    model.add(Dropout(dropout))
    model.add(Dense(1, activation=activation_function_output))
    print('Compiling the CNN-Bidirectional(LSTM) model with two layers of dropout...')
    model.compile(loss=loss, optimizer=optimizer, metrics=['accuracy'])
    return model


def cnn_bilstm_dropout2_gpu(embedding_size, embeddings_initializer, filters, kernel_size,
                            activation_function_conv1d, dropout, max_pool_win, hidden_nodes,
                            activation_function_output, loss, optimizer):
    """This function creates and compile a CNN-Bidirectional(LSTM) model with 2 levels of dropout.
    Note: the model uses the CuDNNLSTM function, which runs only on GPU."""
    print('--- Creating a CNN-Bidirectional(LSTM) model with two layers of dropout (GPU)...')
    model = Sequential()
    model.add(Embedding(NUM_AA + 1, embedding_size, input_length=MAX_PROT_SEQ,
                        embeddings_initializer=embeddings_initializer))
    model.add(Conv1D(filters,
                     kernel_size,
                     padding='valid',
                     activation=activation_function_conv1d,
                     strides=1))
    model.add(MaxPooling1D(pool_size=max_pool_win))
    model.add(Dropout(dropout))
    model.add(Bidirectional(CuDNNLSTM(hidden_nodes)))
    model.add(Dropout(dropout))
    model.add(Dense(1, activation=activation_function_output))
    print('---Compiling the CNN-Bidirectional(LSTM) model with two layers of dropout...')
    model.compile(loss=loss, optimizer=optimizer, metrics=['accuracy'])
    return model
