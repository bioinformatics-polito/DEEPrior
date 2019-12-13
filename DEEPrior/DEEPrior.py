##########################################################################
# DEEPrior: prioritization of gene fusions based on the amino acids
#				  sequence of the resulting proteins
##########################################################################
import click
import pandas as pd
from pathlib import Path
import os
import sys
import re
import keras.backend.tensorflow_backend as ktf

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.create_fusion_list import create_fusion_list
from lib.evaluate_fusions import evaluate_fusions
from lib.model_retraining import model_retraining
from lib.deep_utils import get_session_cpu, get_session_gpu, load_config
from pyfiglet import Figlet
import tensorflow as tf

dir_path = Path(__file__).absolute().parent.as_posix()
results_path = os.path.join(dir_path, "results")
config_path = os.path.join(dir_path, "resources/config.txt")


@click.command()
@click.option('--mode', '-m', help="Tool mode. Only 'retraining' or 'inference' (default) allowed.",
              default='inference')
@click.option('--input_file', '-i',
              help='input file (with path) of gene fusions to prioritize, e.g. '
                   '/home/user/DEEPrior/input_examples/general_out_example.txt.')
@click.option('--fusion_tool', '-f',
              help="Name of the gene fusion detection tool used to obtain the input file. "
                   "Supported tools are: 'DeFuse', 'ChimPIPE', 'EricScript', 'FusionCatcher', 'InFusion', 'JAFFA', "
                   "'SOAPfuse', 'TopHat', 'STAR-Fusion', 'general' (default). To be used only in inference mode.",
              default='general')
@click.option('--version', '-v',
              help="Genome version of input file coordinates. 'grch37' or 'grch38' are allowed.")
@click.option('--model_path', '-mp', help='Path of a retrained model previously obtained with retraining mode. '
                                          'To be specified only in inference mode.',
              default='default')
@click.option('--training', '-t', default='True', help='True to include in the re-training the default training set. '
                                                       'False otherwise. To be used only in retraining mode.')
@click.option('--output', '-o',
              help='Output file (with path and .csv file extension), e.g. /home/user/DEEPrior/results/'
                   'DEEPrior_results.csv. To be used only in inference mode.')
def main(mode, input_file, fusion_tool, version, model_path, training, output):
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    # 1. understand if we are working with GPU or CPu
    IS_GPU = tf.compat.v1.test.is_built_with_cuda()
    if IS_GPU:
        gpu_config = load_config("GPU_PARAMETERS")
        visible_dev_list = gpu_config.get('visible_device_list')

        ktf.set_session(get_session_gpu(float(gpu_config.get('per_process_gpu_memory_fraction')),
                                        bool(gpu_config.get('allow_growth')),
                                        visible_dev_list))
        trained_model_path = os.path.join(dir_path, "resources/HN32_OPTrmsprop_DR2_03_GPU_model.hdf5")
    else:
        cpu_config = load_config("CPU_PARAMETERS")
        ktf.set_session(get_session_cpu(int(cpu_config.get('intra_op_parallelism_threads')),
                                        int(cpu_config.get('inter_op_parallelism_threads'))))
        trained_model_path = os.path.join(dir_path, "resources/HN32_OPTrmsprop_DR2_03_CPU_model.hdf5")

    # 2. print the Logo
    custom_fig = Figlet(font='graffiti')
    print(custom_fig.renderText('DEEPrior'))

    # 3. check the usage mode
    if mode not in ['inference', 'retraining']:
        print("Please check mode parameter. You inserted %s, but allowed options are: inference, retraining" % mode)
        sys.exit()

    # ---------------------------------  INFERENCE MODE  ---------------------------------
    if mode == 'inference':
        # 4a. check if all the parameters required for the inference mode are present: -i, -f, -v, -mp, -o
        if any(n is None for n in [input_file, fusion_tool, version, model_path, output]):
            print("\nSomething went wrong. One or more arguments are missing. Please check that you have entered all"
                  "\nthe arguments correctly (check the file paths, uppercase letters, \
                  \nlowercase letters ...). For more problems take a look at the \
                  \nREADME file.\n\n")
            sys.exit()

        if not os.path.exists(os.path.dirname(input_file)):
            print(
                "Please check input parameter. You inserted an invalid input path. '%s' does not exist." % os.path.dirname(
                    output))
            sys.exit()

        if not os.path.exists(input_file):
            print("Please check input parameter. The '%s' file that does exist." % os.path.basename(input_file))
            sys.exit()

        if fusion_tool not in ['DeFuse', 'STAR-Fusion', 'general', 'ChimPIPE', 'EricScript', 'FusionCatcher',
                               'InFusion', 'JAFFA', 'SOAPfuse', 'TopHat']:
            print("Please check fusion_tool parameter. You inserted '%s', but allowed options are: 'DeFuse', "
                  "'STAR-Fusion', 'ChimPIPE', 'EricScript', 'FusionCatcher', 'InFusion', 'JAFFA', "
                  "'SOAPfuse', 'TopHat', 'general' (default)." % fusion_tool)
            sys.exit()

        if version not in ['grch37', 'grch38']:
            print(
                "Please check version parameter. You inserted '%s', but allowed options are: 'grch37' or 'grch38'" % version)
            sys.exit()

        if model_path == 'default':
            model_path = trained_model_path

        elif not os.path.exists(os.path.dirname(model_path)):
            print(
                "Please check model_path parameter. You inserted an invalid model path. '%s' does not exist."
                % os.path.dirname(model_path))
            sys.exit()

        elif not os.path.exists(model_path):
            print("Please check model_path parameter. The '%s' file that does exist." % os.path.basename(model_path))
            sys.exit()

        elif not os.path.basename(model_path).endswith(".hdf5"):
            print("Please check model path parameter. You inserted an invalid model file. Remember that the model file "
                  "must have a hdf5 extension.")
            sys.exit()

        if not os.path.exists(os.path.dirname(output)):
            print("Please check output parameter. You inserted an invalid output path. '%s' does not exist."
                  % os.path.dirname(output))
            sys.exit()
        if not os.path.basename(output).endswith(".csv"):
            print("Please check output parameter. You inserted an invalid output file. Remember that the output file "
                  "must have a csv extension.")
            sys.exit()

        if os.path.exists(output):
            print("Please check output parameter. You inserted an already existing output file: " + output)
            sys.exit()

        # 5a. if every parameter is correct, open the input file
        try:
            if fusion_tool == 'JAFFA':
                file = pd.read_csv(input_file, sep=',')

            elif fusion_tool == 'TopHat':
                file = pd.read_csv(input_file, sep='\t', header=None)

            else:
                file = pd.read_csv(input_file, sep='\t')
        except IOError:
            print("Problems in opening the input file. Check if it exists and the spelling is correct.")
            sys.exit()

        # 6a. obtain the list of objects with all info about each gene fusion
        list_fusions = create_fusion_list(file, fusion_tool, version, mode)
        # 7a. evaluate the gene fusions
        evaluate_fusions(list_fusions=list_fusions, trained_model_path=model_path, output_file_path=output)

    # ---------------------------------  RETRAINING MODE  ---------------------------------
    elif mode == 'retraining':
        # 4a. check if all the required parameters for the retraining mode are present: -i, -v, -t
        if any(n is None for n in [input_file, version, training]):
            print("\nSomething went wrong. One or more arguments are missing. Please check that you have entered all \
                   \nthe arguments correctly (check the file paths, uppercase letters, "
                  "\nlowercase letters ...). For more problems take a look at the "
                  "\nREADME file.\n\n")
            sys.exit()

        if not os.path.exists(os.path.dirname(input_file)):
            print("Please check input parameter. You inserted an invalid input path. '%s' does not exist."
                  % os.path.dirname(output))
            sys.exit()

        if not os.path.exists(input_file):
            print("Please check input parameter. The '%s' file that does exist." % os.path.basename(input_file))
            sys.exit()

        if version not in ['grch37', 'grch38']:
            print("Please check version parameter. You inserted '%s', but allowed options are: 'grch37' or 'grch38'."
                  % version)
            sys.exit()

        if training not in ['True', 'False']:
            print(
                "Please check the training parameters. You inserted '%s', but allowed options are: 'True' (default) or 'False'."
                % training)
            sys.exit()
        elif training == 'True':
            training_flag = True
        elif training == 'False':
            training_flag = False
        # 5a. if every parameter is correct, open the input file
        try:
            file = pd.read_csv(input_file, sep='\t')
        except IOError:
            print("Problems in opening the input file. Check if it exists and the spelling is correct.")
            sys.exit()

        # 6a. obtain list of objects with all info about each gene fusion
        list_fusions = create_fusion_list(file, fusion_tool, version, mode)
        print("\n\nDatabase succesfully built, evaluating fusions with Deep Learning model...")

        # 7a. retrain the model with the new training set
        model_retraining(fusions_list=list_fusions, add_default=training_flag, is_gpu=IS_GPU)

    return


if __name__ == '__main__':
    main()
