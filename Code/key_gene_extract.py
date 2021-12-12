import miRNET
import json

Targets = miRNET.Targets(path_to_miRTarBase='./baseData/hsa_miRTarBase.csv')

miR_key_dict = dict()
tis_gene_set = miRNET.tissue_selector(ans=0, tissue_id=23)

for miR in list(Targets.miR_dict.keys()):
#    if miR[-1] == 'p':
#        miR = miR[:-3] + '-'
    miR_targets = Targets.get_targets(miR, False)
    if miR_targets == 1:
        continue
    MirNet = miRNET.MainNet()  # Load String db and create gene-gene interaction network
    MirNet.get_LCC()  # get the largest connected component from the network
    check = 0
    check += MirNet.select_nodes(miR_targets)  # select the part of LCC containing only the miRNA target genes
    check += MirNet.select_nodes(tis_gene_set)  # select the part of LCC containing only the tissue target genes
    if check == 1:
        continue
    kne = miRNET.KeyNodesExtractor(MirNet)
    kne_ex = kne.extraction()
    if kne_ex != 1:
        key_genes = list(kne_ex)
    else:
        continue
    miR_key_dict[miR] = {'targets': len(miR_targets),
                         'LCC': len(MirNet.LCC.nodes()),
                         'key_genes': key_genes}

with open('./addData/miR_key_dict.json', 'w') as outfile:
    json.dump(miR_key_dict, outfile)
