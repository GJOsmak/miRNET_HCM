import pandas as pd
import numpy as np
import scipy as sp
import miRNET, miRNET_enrichment
#Targets = miRNET.Targets(path_to_miRTarBase='./baseData/hsa_miRTarBase.csv')
import warnings
import collections
from matplotlib import pyplot as plt
import random
from scipy.stats import chi2_contingency
import json
import networkx as nx

with open('./addData/miR_key_dict.json') as json_file:
    miR_key_dict = json.load(json_file)

dict_for_mc = dict()

for miR in miR_key_dict.keys():
    if 1351 > miR_key_dict[miR]['targets'] > 51:
        if miR_key_dict[miR]['LCC'] > 5:
            if len(miR_key_dict[miR]['key_genes']) > 5:
                dict_for_mc[miR] = miR_key_dict[miR]

dict_miR_to_paths = dict()
for miR, genes in dict_for_mc.items():
    print(miR)
    enrich_res = miRNET_enrichment.reactome_enrichment(genes['key_genes'], species='Homo sapiens')
    enrich_res = miRNET_enrichment.reac_pars(enrich_res)
    dict_miR_to_paths[miR] = list(enrich_res.react_dict.keys())
with open('./addData/miR_path_dict.json', 'w') as outfile:
    json.dump(dict_miR_to_paths, outfile)