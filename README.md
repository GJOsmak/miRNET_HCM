# miRNET_HCM

Step-by-step analysis of the function of the HCM-associated microRNAs expressed in heart tissue and described in the [review](https://www.internationaljournalofcardiology.com/article/S0167-5273(20)34078-X/fulltext)

## 1. Data preparation

Execute key_gene_extract.py for creating the dictionary {miRNA: {targets=int(), LCC=int(), key_genes=list()}}, where targets - is a total amount of the miRNA targets in miRTraBase; LCC - is the size of the Largest Connected Component of a miRNA target-gene  network; key_genes - list of key genes selected from LCC using [miRNET algorithm](https://pubmed.ncbi.nlm.nih.gov/33352947/).

This dictionary will be used for Monte Carlo tests in step #4.

## 2. Bigraph "miRNA" -> "Reactome gene set" constaraction

For this use bigraph_construction.ipynb

## 3. Bigraph "miRNA" -> "Reactome gene set" constaraction

For the down-regulated miRNAs the bigraph was disconnected and we analyse them separetly, one by one.
This analysis is performed in down_reg_miRNA_one_by_one.ipynb

## 4. Monte Carlo tests

The Monte Carlo tests of significant of the results from #3 and #4 is performed in MC.ipynb

## 5. To study the regulatory effects of up-miRNA on the pathways from bigraph #2

Use up_reg_miR_by_paths.ipynb