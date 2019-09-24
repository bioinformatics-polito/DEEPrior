from keras.callbacks import EarlyStopping, ModelCheckpoint
import sys,os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from pathlib import Path
from lib.deep_utils import create_val_train_sets, load_tokenizer, load_config, cnn_bilstm_dropout2_cpu, \
    cnn_bilstm_dropout2_gpu

ONCO = 1
NON_ONCO = 0


def model_retraining(fusions_list, add_default, is_gpu):
    dir_path = Path(__file__).absolute().parent.parent.as_posix()
    output_path = os.path.join(dir_path, 'results')
    best_model_file_path = os.path.join(output_path, "model.hdf5")

    # 1. load the tokenizer
    tokenizer = load_tokenizer()

    # 2. read the callback and the model parameters from the config file
    config_callback = load_config("CALLBACK_PARAMETERS")
    print("Retrieving the model parameter from the configuration file..")
    if config_callback is None:
        print("Error in the configuration file. No section CALLBACK_PARAMETERS")

    config_model = load_config("MODEL_PARAMETERS")
    if config_model is None:
        print("Error in the configuration file. No section MODEL_PARAMETERS")

    # 3. create the validation and the training set from the fusion list
    print("Creating a validation and a training set..")
    x_train, y_train, x_val, y_val = create_val_train_sets(config_model.get("validation_set_percentage"),
                                                           fusions_list,
                                                           tokenizer,
                                                           add_default)
    # 4. set callbacks parameters
    callbacks = [EarlyStopping(monitor=config_callback.get("earlystopping_monitor"),
                               patience=int(config_callback.get("patience")),
                               min_delta=float(config_callback.get("min_delta")),
                               verbose=1),
                 ModelCheckpoint(filepath=best_model_file_path,
                                 monitor=config_callback.get("modelcheckpoint_monitor"),
                                 save_best_only=config_callback.get("save_best_only"),
                                 mode=config_callback.get("mode"),
                                 period=int(config_callback.get("period")))]

    # 5. create the model with config params
    if is_gpu:
        model = cnn_bilstm_dropout2_gpu(embedding_size=int(config_model.get("embedding_size")),
                                        embeddings_initializer=config_model.get("embeddings_initializer"),
                                        filters=int(config_model.get("filters")),
                                        kernel_size=int(config_model.get("kernel_size")),
                                        activation_function_conv1d=config_model.get("activation_function_conv1d"),
                                        dropout=float(config_model.get("dropout")),
                                        max_pool_win=int(config_model.get("max_pooling_window")),
                                        hidden_nodes=int(config_model.get("hidden_nodes")),
                                        activation_function_output=config_model.get("activation_function_output"),
                                        loss=config_model.get("loss"),
                                        optimizer=config_model.get("optimizer"))
    else:
        model = cnn_bilstm_dropout2_cpu(embedding_size=int(config_model.get("embedding_size")),
                                        embeddings_initializer=config_model.get("embeddings_initializer"),
                                        filters=int(config_model.get("filters")),
                                        kernel_size=int(config_model.get("kernel_size")),
                                        activation_function_conv1d=config_model.get("activation_function_conv1d"),
                                        dropout=float(config_model.get("dropout")),
                                        max_pool_win=int(config_model.get("max_pooling_window")),
                                        hidden_nodes=int(config_model.get("hidden_nodes")),
                                        activation_function_output=config_model.get("activation_function_output"),
                                        loss=config_model.get("loss"),
                                        optimizer=config_model.get("optimizer"))
    # 6. fit the model
    print('Fitting the model...')
    model.fit(x_train, y_train, batch_size=int(config_model.get("batch_size")),
              epochs=int(config_model.get("epochs")), validation_data=[x_val, y_val],
              callbacks=callbacks)

    print("The trained the model is saved at: "+best_model_file_path)

