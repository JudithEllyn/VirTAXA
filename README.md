# VirTAXA

VirTAXA is a robust tool to conduct taxonomic classification of RNA viruses. To use VirTAXA, you only need to input your contigs and their proteins to the program.

## An easy way to install

- We recommend you to install all the package with Anaconda
  
  - After downloading this respository, you can use anaconda to install the virtaxa.yaml. This will install all packages you need.
  
  ```
  conda env create -f virtaxa.yaml -n virtaxa
  conda activate virtaxa
  ```
  
  - Noting that sometimes you need to install [TreeShrink](https://github.com/uym2/TreeShrink) manually to make sure it works (use `run_treeshrink.py -h` to test).
- If you do not want to use the .yaml, you can also install the requirements below independently by yourself.

## Requirements

- Required
  
  - Python3
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

- There are 3 required parameters for the program: 1. --query/-q protein file of the query sequences. 2. --db/-d the database (default **db_vir**  we build) 3. --output/-o the folder saving the prediction results.
  <br/>
- Here we present an example to show how to run VirTAXA. We provide a test file named “vir_test.fasta" and its corresponding protein file "vir_test.faa" in the Github folder. The only command that you need to run is as below. However, before using the database, do not forget to download the db_vir as the database and unzip the **hmm_train.hmm** and **train.csv** files.
  
  ```python Virtaxa.py --db db_vir --query vir_test.faa --output $output_folder$ ```
  <br/>
- The output file is **final_result.csv** in the result folder. There are three column in this csv file: "Accession, pred_genus, step". The “Step” column provides the user with the steps leading to this prediction so that the user can decide whether to use all or part of the prediction results. **high_level.csv** is the result on higher level for the sequences failed to classified at the genus level.
- The **aln_file** folder occupies a large amount of space . If you do not use the `--fast` option, there is no need to download this folder.
  <br/>
- There are 2 optional parameters.
  
  - --new. This parameter is used to predict possible new genus-level clusters using the phylogenetic tree analysis. This functionality requires additional processing time. Users can leverage this function by adding this parameter in the command. The result file  **new_cluster.json** provides the putative new genus-level clusters. Note that these new clusters may not necessarily correspond to strict genus-level units in the stantard taxonomy. VirTAXA just provides them for a reference.
  - --fast. This is a new version that can speed up the multiple sequence alignment (MSA) processing by using the pre-built alignments. Users can leverage this function by adding this parameter in the command.  Need to download the **aln_file** folder in the reference database before using this function.

## Notice

- If you want to use VirTAXA, you need to take care of these things:

1. Make sure all your contigs are virus contigs, or they can create erroneous results.
2. When using your own query data, nucleotide sequences (.fasta) and protein sequences (.faa) are both required and make sure they have the same prefix, for example, vir_test.
