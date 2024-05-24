# VirTAXA

VirTAXA is a robust tool to conduct taxonomic classification of RNA viruses. To use VirTAXA, you only need to input your contigs and their proteins to the program.


## An easy way to install
- We recommend you to install all the package with Anaconda

    - After cloning this respository, you can use anaconda to install the environment.yaml. This will install all packages you need with gpu mode (make sure you have installed cuda on your system).

    ```
    conda env create -f environment.yaml -n virtaxa
    conda activate virtaxa
    ```
- If not, you can also install the requirements below independently by yourself. 

## Requirements
- Required
    -  Python3
    - Diamond
    - BLAST
    - HMMER3
    - FastTree
    - TreeShrink
    - MAFFT

- Python3 Dependencies
    - Pandas
    - Numpy
    - Biopython
    - tqdm

## Usage (example)
- Here we present an example to show how to run VirTAXA. We support a file named “test.fasta" and its corresponding protein file "test.faa" in the Github folder for testing the software. The only command that you need to run is as below. However, before using the database, do not forget to download the db_vir as the database and extract the hmm_train.hmm and train.fasta files in it.

    ```python Virtaxa.py --db db_vir --query vir_test.faa --output $output_folder$ ```

- There are 3 parameters for the program: 1. --query/-q protein file of your query sequences. 2. --db/-d the database (default db_vir which we build us) 3. --output/-o the folder you want to save the prediction results.

- The output file is final_result.csv in the result folder. There are three column in this csv file: "Accession, pred_genus, step". "step" column provides the users with which step leads to this prediction, allowing users to decide whether to include all or part of the classification outputs.


## Notice
- If you want to use VirTAXA, you need to take care of these things:
1. Make sure all your contigs are virus contigs, or they can create erroneous results.
2. When you use your own query data, test.fasta and test.faa are both required and make sure they have the same prefix.


