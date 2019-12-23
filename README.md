# DEEPrior

DEEPrior is an inherently flexible deep learning tool that predicts the probability
of a gene fusion being a driver of an oncogenic process, by directly exploiting the amino acid sequence of the fused protein and it is able to prioritize gene fusions from different tumors. Unlike state of the art tools, it also support easy retraining and re-adaptation of the model. 

To implement these concepts, the tool exploits two modes: inference (to prioritize gene fusions) and retraining (to create a new deep learning model).

DEEPrior is implemented in Python 3.7 with minimal additional libraries, and it is available both for CPU and GPU.

In the following you will find:

1) *Getting Started*: obtain a working copy of DEEPrior
2) *Usage*: how to use DEEPrior with examples
3) *Files*: input and output files for inference and retraining mode
4) *Datasets*: description of datasets used to train and test DEEPrior


## 1. Getting Started

### 1.1 Prerequisites
DEEPrior is developed in Python 3.7 with minimal libraries required. To run DEEPrior we strongly suggest you to create a clean virtual environment in order to avoid conflicts with other projects. If you are an expert with virtual environments, all you need is to install the libraries listed in the requirements files, clone this repository and jump directly to **Test if everything is ok**. Otherwise no problem, follow the **Installing** section, the installation is very simple!

#### 1.1.1 Prerequisites CPU

- xlrd 1.2.0
- pandas 0.24.2
- biopython 1.73
- matplotlib 3.1.0
- numpy 1.16.4
- requests 2.22.0
- scikit-learn 0.21.2
- scipy 1.2.1
- click 7.0
- configparser 3.7.4
- pyfiglet 0.8.post1
- keras 2.2.4
- tensorflow 1.13.1

The prerequisites are listed in the requirements_CPU.txt file. 

#### 1.1.2 Prerequisites GPU
We assume that **cuda 10.0** is installed on your system.

- xlrd 1.2.0
- pandas 0.24.2
- biopython 1.73
- matplotlib 3.1.0
- numpy 1.16.4
- requests 2.22.0
- scikit-learn 0.21.2
- scipy 1.2.1
- click 7.0
- configparser 3.7.4
- pyfiglet 0.8.post1
- keras 2.2.4
- tensorflow 1.13.1
- tensorflow-estimator 1.13.0 
- tensorflow-gpu 1.13.1

The prerequisites are listed in the requirements_GPU.txt file. 

### 1.2 Installing
First of all, check if you have pip and the virtual environments packages for Python3. If pip and/or virtualenv are not installed in your system, follow the instructions reported [here](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/).

Now copy and paste these instrucions to create and activate a DEEPrior virtual environment called *DEEPenv*:
```
python3 -m venv DEEPenv
source DEEPenv/bin/activate
```
Then move in the DEEPenv forlder, **clone this repository** and install all the required packages. 
Use requirements_CPU.txt if you are installing the CPU version of DEEPrior, requirements_GPU.txt otherwise. 
```
cd DEEPenv
git clone https://github.com/bioinformatics-polito/DEEPrior.git
cd DEEPrior
pip3 install -r requirements_CPU.txt
```

### 1.3 Test if everything is ok
Once you have followed the previous steps, **move into DEEPenv folder** and test the tool with the following commands:

```
cd ../              # move to DEEPenv folder. You can use the global path e.g. cd /home/user/DEEPenv
source bin/activate # command to activate virtual environment if you followed our installation guide
cd DEEPrior/DEEPrior
mkdir results
python DEEPrior.py -i input_examples/general_out_example.txt -f general -v grch37 -o results/DEEPrior_results.csv
```

If everything worked correctly, after a few seconds you will find DEEPrior_results.csv file in the DEEPenv/DEEPrior/DEEPrior/results folder.
It's done, you are now ready to use DEEPrior!


## 2. Usage

DEEPrior has two main modes: inference and retraining. The mode is specified through the -m --mode paramter. By default, it is equal to 'inference'.

  **-h, --help**               Shows this help message and exit
  **-m --mode**                Tool mode: 'retraining' or 'inference' (default) allowed.


### 2.1 Inference mode

The usage of the tool in inference mode is the following one:
```
python DEEPrior.py [-h] [-m inference] [-i INPUT] [-f FUSION_TOOL] [-v VERSION] [-mp MODEL_PATH] [-o OUTPUT]
```

  **-i INPUT, --input INPUT**
                        input file (with path) of gene fusions to prioritize,
                        e.g. /home/user/DEEPenv/DEEPrior/DEEPrior/input_examples/general_out_example.txt

  **-f FUSION_TOOL, --fusion_tool FUSION_TOOL**
                        Name of the gene fusion detection tool used to obtain
                        the input file. Supported tools are: 'DeFuse', 'STAR-Fusion', 'ChimPIPE', 'EricScript', 'FusionCatcher',
                        'InFusion', 'JAFFA', 'SOAPfuse', 'TopHat', 'general' (default)

  **-v VERSION, --version VERSION**
                        Genome version of input file coordinates. 'grch37' or
                        'grch38' allowed

  **-mp MODEL_PATH, --model_path MODEL_PATH**
                        Path that points to the model to be used for the inference, to be specified only if a new model has been obtained with *retraining* mode
                        If 'default' or not specified, the native model is used.

  **-o OUTPUT, --output OUTPUT**
                        Name (with path) of the output file, e.g.
                        /home/user/DEEPenv/DEEPrior/DEEPrior/results/DEEPrior_results.csv 
                        The output file extension must be .csv

Example:
```
python DEEPrior.py -i /home/user/DEEPenv/DEEPprior/DEEPrior/input_examples/general_out_example.txt -f general -v grch37 -o /home/user/DEEPenv/DEEPrior/DEEPrior/results/DEEPrior_results.csv
```
Please, have a look to **Input files** and **Output file** for details about the files.


### 2.2 Reatraining mode
The usage of the tool in inference mode is the following one:
```
python DEEPrior.py [-h] [-m retraining] [-i INPUT] [-v VERSION] [-t TRAINING_FLAG]
```

  **-i INPUT, --input INPUT**
                        input file (with path) of gene fusions to prioritize,
                        e.g. /home/user/DEEPenv/DEEPrior/DEEPrior/input_examples/re_train_example.csv

  **-v VERSION, --version VERSION**
                        Genome version of input file coordinates. 'grch37' or
                        'grch38' allowed

  **-t TRAINING_FLAG --training TRAINING_FLAG**
                        True to include in the re-training process the default training set. 
                        False otherwise.


## 3 Files

### 3.1 Inference mode
Inference mode is the default mode. The following are the input and output files.

#### 3.1.1 Inference mode input file
DEEPrior is intended to be run after gene fusion detection tools in order to prioritize the output and focus on gene fusions with a higher probability to be involved in oncogenic processes.
DEEPrior natively support the output of DeFuse, STAR-Fusion, ChimPIPE, EricScript, FusionCatcher, InFusion, JAFFA, SOAPfuse, TopHat, however any gene fusion can be processed providing in a **tab separated file** the genomic coordinates of the breakpoints (see *general_out_example.txt* file in *input_example* folder).

An example of the *general* format is the following:

| chr5p | coord5p | chr3p | coor3p   |
|-------|---------|-------|----------|
| chr7  | 1000000 | chr4  | 1000000  |
| chr9  | 2555965 | chr6  | 56444888 |


The first two columns refer to chromosome number and breakpoint coordinate of 5p gene, while third and fourth columns refer to 3p gene. Coordinates can be entered in genome version *grch37* or *grch38*.

In *input_examples* folder you can find examples of all the allowed input files (general and gene fusion detection tools output).

#### 3.1.2 Inference mode output file
The output file contains the following information

- **fusion_pair:** name of the gene fusion with common gene names 
- **oncogenic probability value:** oncogenic probability value reported by the tool. It is a number between 0 and 1. Closer is the number to 1, higher is the probability to be oncogenic
- **version:** grch37 or grch38 depending on the genome version parameter defined during the running of DEEPrior. Remember that hg19 is equivalent to grch37 and hg38 is equivalent to grch38
- **chr5p:** chromosome number of 5p gene
- **coord5p:** breakpoint coordinate of 5p gene on chromosome 5p (1-based coordinate system)
- **5p strand:** strand of 5p gene
- **5p common name:** common name of 5p gene
- **5p ensg:** ENSEMBL gene identifier of 5p gene
- **5p gene functionality:** functionality of 5p gene (e.g. proteing coding or not)
- **5p gene description:** additional information about 5p gene provided by ENSEMBL, usually a description of the biological process in which the gene is involved
- **chr3p:** chromosome number of 3p gene
- **coord3p:** breakpoint coordinate of 3p gene on chromosome 3p (1-based coordinate system)
- **3p strand:** strand of 3p gene
- **3p common name:** common name of 3p gene
- **3p ensg:** ENSEMBL gene identifier of 3p gene
- **3p gene functionality:** functionality of 3p gene (e.g. proteing coding or not)
- **3p gene description:** additional information about 3p gene provided by ENSEMBL, usually a description of the biological process in which the gene is involved
- **MainProteinLength:** length of the fused protein
- **TruncatedProtein:** Yes if the fused protein is truncated (an early stop codon occurs in the protein). No otherwise.
- **5p_gene_complete:** Yes if 5p gene is complete in the fusion (stop codon in upstream gene is present in the protein). No otherwise.
- **3p_gene_complete:** Yes if 3p gene is complete in the fusion (start codon in downstream gene is present in the protein). No otherwise.
- **main protein:** the protein with no skipped exons


### 3.2 Retraining mode
Although the retraining mode is not the main one, DEEPrior allows you to retrain the deep learning model if new validated gene fusions are available.

#### 3.2.1 Retraining mode input file
In this case, the input file is a **tab separated file** and contains validated gene fusions to be included in the retraining of the model for which the label (*oncogenic* or *not oncogenic*) is known.
The file is similar to the one reported in 3.1.1 and in addition it contains the **Label** column which indicates the class to which that gene fusion belongs. 0 means not oncogenic and 1 oncogenic. An example of the file is provided below and also in the *input_examples* folder (*re_train_example.csv* file):

| chr5p | coord5p   | chr3p | coor3p    | label |
|-------|-----------|-------|-----------|-------|
| 6     | 31637695  | 22    | 42486683  | 1     |
| 2     | 85132848  | 2     | 37456103  | 1     |
| 19    | 14676464  | 16    | 10862958  | 1     |
| 21    | 47542847  | 16    | 30581751  | 0     |
| 5     | 134261417 | 2     | 219144833 | 0     |
| 5     | 134262389 | 2     | 1499769   | 0     |


#### 3.2.2 Retraining mode output file
The retraining mode output consists of a **.hdf5 file** containing the weights and the architecture of the new trained model. This model can then be used to make the gene fusions inference instead of the default deep learning model.

## 4. Datasets
Together with the DEEPrior tool, we provide to the scientific community the datasets used to train and test DEEPrior. The datasets (*training.csv*, *test_set_1.csv*, *test_set_2.csv*) are available in *DEEPrior/data* folder in the format described in section 2.2. Moreover, there is a *Label* column identifying which class the gene fusion belongs to (0 for *NotOnco* and 1 for *Onco*).

### 4.1 Training set
This set consists of 786 fusion pairs (777 genes overall) and 2118 sequences, respectively 1059 Onco and 1059 NotOnco, obtained from COSMIC (Catalog of Somatic Mutations in Cancer) and Babicenau et al. work *Recurrent chimeric fusion rnas in non-cancer tissues and cells*.

### 4.2 Test set 1
This set is composed of a total of 142 fusion pairs and 156 gene fusions, of which 122 Onco and 34 NotOnco. The sequences associated with Onco gene fusions were extracted from the ChimerDB2.0 database, while NotOnco were constituted by the false positives reported by TopHat-Fusion and STAR-Fusion on Illumina BodyMap 2.0 samples.

### 4.3 Test set 2
This set is composed of 2595 fusion pairs and 2623 gene fusions, all belonging to the Onco category. This dataset was built starting from the work of Gao et al. *Driver fusions and their implications in the development and treatment of human cancers*.

## Authors

* **Marta Lovino** - contact marta.lovino@polito.it
* **Maria Serena Ciaburri** 
* **Gianvito Urgese** 
* **Enrico Macii** 
* **Santa Di Cataldo** 
* **Elisa Ficarra** - contact elisa.ficarra@polito.it


## License

This project is licensed under the AGPL v3 License - see the [LICENSE.rst](LICENSE.rst) file for details

## Acknowledgments

We thank Gao et al., 2018 (*Driver fusions and their implications in the development and treatment of human cancers*) for providing WGS validated data and Wen-Wei Liang for illustrating the WGS validation process.

Furthermore, a special thanks to all those who will use DEEPrior and will provide suggestions for improvement.

