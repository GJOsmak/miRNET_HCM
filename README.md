# miRNET_HCM

Step-by-step analysis of the function of the HCM-associated microRNAs expressed in heart tissue

## 1. Data preparation

Execute `key_gene_extract.py` for creating the dictionary {miRNA: {targets=int(), LCC=int(), key_genes=list()}}, where targets - is a total amount of the miRNA targets in miRTraBase; LCC - is the size of the Largest Connected Component of a miRNA target-gene  network; key_genes - list of key genes selected from LCC using [miRNET algorithm](https://pubmed.ncbi.nlm.nih.gov/33352947/).

This dictionary will be used for Monte Carlo tests in step #4.

Execute `miR_to_paths_creator.py` for creating the dictionary {miRNA: list(Reactome pathways overrepresented in key genes)}.

This dictionary will be used for `pathway_test.ipynb` in step #5.

## 2. Bigraph "miRNA" -> "Reactome gene set" constaraction

For this use `bigraph_construction.ipynb`

## 4. Monte Carlo tests

The Monte Carlo tests of significant of the results is performed in `MC.ipynb`

## 5. Estimation of the probability of obtaining the signaling pathways

Use `pathway_test.ipynb`

## 6. To study the regulatory effects of up-miRNA on the pathways from bigraph #2

Use `up_reg_miR_by_paths.ipynb`

## 7. Analysis of the consistency of the algorithm used to select key target genes 

Use `Node_test.ipynb`