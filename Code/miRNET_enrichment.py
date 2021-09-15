from reactome2py import analysis
import requests  # HTTP Client for Python
import json  # Standard JSON library
from py2cytoscape import util as cy
import networkx as nx
import pandas as pd
import io
import itertools as itt
from scipy.spatial.distance import pdist
import scipy.cluster
import matplotlib.pyplot as plt
from colour import Color

def reactome_enrichment(gene_set, species='Homo sapiens'):
    target_set = ",".join(gene_set)

    result = analysis.identifiers(ids=target_set)
    token = result['summary']['token']

    token_result = analysis.token(token, species=species, page_size='-1', page='-1', sort_by='ENTITIES_FDR',
                                  order='ASC', resource='TOTAL', p_value='0.05', include_disease=False,
                                  min_entities=2, max_entities=None)

    return token_result


def reac_pars(token_result):

    class ReactomeRes:
        pass

    ReactomeRes.token = token_result['summary']['token']
    ReactomeRes.max_p_value = 0
    ReactomeRes.min_p_value = 99999
    ReactomeRes.max_ratio = 0
    ReactomeRes.min_ratio = 99999
    ReactomeRes.min_found = 99999
    ReactomeRes.max_found = 0
    ReactomeRes.react_dict = dict()

    for item in token_result['pathways']:

        if item['entities']['found'] < 2:
            continue

        ReactomeRes.react_dict[item['stId']] = {'name': item['name'],
                                      'fdr': item['entities']['fdr'],
                                      'p_value': item['entities']['pValue'],
                                      'found': item['entities']['found'],
                                      'ration': item['entities']['ratio']}

        #   for Networks

        if item['entities']['found'] < ReactomeRes.min_found:
            ReactomeRes.min_found = item['entities']['found']
        if item['entities']['found'] > ReactomeRes.max_found:
            ReactomeRes.max_found = item['entities']['found']

        if item['entities']['ratio'] < ReactomeRes.min_ratio:
            ReactomeRes.min_ratio = item['entities']['ratio']
        if item['entities']['ratio'] > ReactomeRes.max_ratio:
            ReactomeRes.max_ratio = item['entities']['ratio']

        if item['entities']['pValue'] < ReactomeRes.min_p_value:
            ReactomeRes.min_p_value = item['entities']['ratio']
        if item['entities']['pValue'] > ReactomeRes.max_p_value:
            ReactomeRes.max_p_value = item['entities']['ratio']

    return ReactomeRes


def get_net(ReactomeRes):

    reactome_linkege = pd.read_csv('./baseData/ReactomePathwaysRelation.txt', sep='\t', header=None)
    reactome_linkege.columns = ['Source', 'Target']

    G = nx.from_pandas_edgelist(reactome_linkege, source='Source', target='Target')
    G_enrich = G.subgraph(ReactomeRes.react_dict.keys())

    # removing network components that contain less than two nodes

    G_enrich = nx.Graph(G_enrich)

    for CC in list(nx.connected_components(G_enrich)):
        if len(CC) < 3:
            G_enrich.remove_nodes_from(CC)

    # add Statistics to Network

    for nods, dct in G_enrich.nodes(data=True):
        dct['label'] = ReactomeRes.react_dict[nods]['name']
        dct['fdr'] = ReactomeRes.react_dict[nods]['fdr']
        dct['found'] = ReactomeRes.react_dict[nods]['found']
        dct['ration'] = ReactomeRes.react_dict[nods]['ration']
        dct['p_value'] = ReactomeRes.react_dict[nods]['p_value']

    return G_enrich


def dendro_reactome_to_pandas(ReactomeRes, G, species='Homo sapiens'):

    url = 'https://reactome.org/AnalysisService/download/' + ReactomeRes.token + \
          '/pathways/TOTAL/result.csv'
    res = requests.get(url).content
    reactome_df = pd.read_csv(io.StringIO(res.decode('utf-8')))

    reactome_df = reactome_df[reactome_df['Entities pValue'] < 0.05]
    reactome_df = reactome_df[reactome_df['#Entities found'] >= 2]
    reactome_df = reactome_df[reactome_df['Species name'] == species]
    reactome_df = reactome_df[reactome_df['Pathway identifier'].isin(G.nodes)]

    return reactome_df

"""
def draw_net_to_cytoscape(G, ReactomeRes):

    PORT_NUMBER = 1234
    IP = 'localhost'
    BASE = 'http://' + IP + ':' + str(PORT_NUMBER) + '/v1/'

    # requests.delete(BASE + 'session')  # Delete all networks in current session

    cytoscape_network = cy.from_networkx(G)
    cytoscape_network['data']['name'] = 'First_Enrich_Net'
    res1 = requests.post(BASE + 'networks', data=json.dumps(cytoscape_network))
    res1_dict = res1.json()
    new_suid = res1_dict['networkSUID']
    requests.get(BASE + 'apply/layouts/force-directed/' + str(new_suid))

    # load and apply style

    res = requests.get(BASE + 'styles/PathwayEnrichStyles')
    if res.status_code != 200:

        with open('./options/cytoscape_styles/PathwayEnrichStyles.json') as json_file:
            directed_styles = json.load(json_file)

        for mapings in range(0, len(directed_styles['mappings'])):
            if directed_styles['mappings'][mapings]['visualProperty'] == 'NODE_FILL_COLOR':
                directed_styles['mappings'][mapings]['points'][0]['value'] = ReactomeRes.min_p_value
                directed_styles['mappings'][mapings]['points'][2]['value'] = 0.05
            if directed_styles['mappings'][mapings]['visualProperty'] == 'NODE_TRANSPARENCY':
                directed_styles['mappings'][mapings]['points'][0]['value'] = ReactomeRes.min_found
                directed_styles['mappings'][mapings]['points'][1]['value'] = ReactomeRes.max_found
            if directed_styles['mappings'][mapings]['visualProperty'] == 'NODE_LABEL_TRANSPARENCY':
                directed_styles['mappings'][mapings]['points'][0]['value'] = ReactomeRes.min_found
                directed_styles['mappings'][mapings]['points'][1]['value'] = ReactomeRes.max_found
            if directed_styles['mappings'][mapings]['visualProperty'] == 'NODE_BORDER_TRANSPARENCY':
                directed_styles['mappings'][mapings]['points'][0]['value'] = ReactomeRes.min_found
                directed_styles['mappings'][mapings]['points'][1]['value'] = ReactomeRes.max_found
            if directed_styles['mappings'][mapings]['visualProperty'] == 'NODE_SIZE':
                directed_styles['mappings'][mapings]['points'][0]['value'] = ReactomeRes.min_ratio
                directed_styles['mappings'][mapings]['points'][1]['value'] = ReactomeRes.max_ratio

        # Create new Visual Style
        res = requests.post(BASE + "styles", data=json.dumps(directed_styles))
        # new_style_name = res.json()['title']

        # Apply it to current network

    requests.get(
        BASE + 'apply/styles/' + 'PathwayEnrichStyles' + '/' + str(new_suid))  # !Это говно почему-то не работает

    return res
"""

# Dendrograms


def create_similarity_matrix(gene_sets):
        """Create a similarity matrix for a given pathway-geneset dataset.
        :param dict gene_sets: pathway gene set dictionary
        :rtype: pandas.DataFrame
        :returns: similarity matrix
        """
        index = sorted(gene_sets.keys())
        similarity_dataframe = pd.DataFrame(0.0, index=index, columns=index)

        for pathway_1, pathway_2 in itt.product(index, index):
            intersection = len(gene_sets[pathway_1].intersection(gene_sets[pathway_2]))
            smaller_set = min(len(gene_sets[pathway_1]), len(gene_sets[pathway_2]))

            similarity = float(intersection / smaller_set)  # Formula to calculate similarity

            similarity_dataframe[pathway_1][pathway_2] = similarity

        return similarity_dataframe


def get_dendro(dt, key_nodes_dict,mir_name):

    gene_set_dict = {dt.iloc[row]['Pathway identifier']: set(dt.iloc[row]['Submitted entities found'].split(';')) for
                     row in range(0, dt.shape[0])}
    dtid2name = {dt.iloc[row]['Pathway identifier']: set(dt.iloc[row]['Pathway name'].split(';')) for row in
                 range(0, dt.shape[0])}

    similarity_matrix = create_similarity_matrix(gene_set_dict)
    
    print(similarity_matrix)

    distance_matrix = pdist(similarity_matrix, metric='correlation')
    distance_matrix = pd.Series(distance_matrix).fillna(1)
    
    print(distance_matrix)
    print(distance_matrix.shape)

    # Cluster hierarchicaly using scipy
    clusters = scipy.cluster.hierarchy.linkage(distance_matrix, method='single')
    T = scipy.cluster.hierarchy.to_tree(clusters, rd=False)

    # Create dictionary for labeling nodes by their IDs
    labels = list(similarity_matrix.columns)
    id2name = dict(zip(range(len(labels)), labels))

    def getlabel(x):
        return '{}: ({})'.format(''.join(dtid2name[id2name[x]]), ', '.join(gene_set_dict[id2name[x]]))

    # Draw dendrogram using matplotlib to scipy-dendrogram.pdf
    my_dpi = 10
    plt.figure(figsize=(2400 / my_dpi, 2400 / my_dpi), dpi=my_dpi)
    fig, axes = plt.subplots()
    scipy.cluster.hierarchy.dendrogram(clusters,
                                       leaf_label_func=getlabel,
                                       orientation='right')
    fig.set_figwidth(10)
    fig.set_figheight(16)
    plt.tight_layout()

    """
    set a label color as the normalised cumulative betweenness {sum(betweenness)/len(gene_set)}
    """
    top5centralNods = key_nodes_dict

    ylbls = axes.get_ymajorticklabels()

    # create label to normalised cumulative betweenness dictionary

    label_to_BC = dict()
    for lbl in ylbls:
        local_label = (lbl.get_text().split(': ('))[0]  # get path name
        gene_set = (lbl.get_text().split(': ('))[1][:-1].split(', ')  # get gene set
        sum_BC = sum([top5centralNods[gene] for gene in gene_set]) / len(gene_set)
        label_to_BC[local_label] = sum_BC

    label_to_BC_sort = {k: v for k, v in sorted(label_to_BC.items(), key=lambda item: item[1], reverse=False)}

    # create gradient (grey_to_red hist path
    grey = Color('grey')
    colors = list(grey.range_to(Color("red"), len(label_to_BC_sort)))
    for label, tmp_color in zip(label_to_BC_sort.keys(), colors):
        color = str(tmp_color)
        if len(color) < 7 and color[0] == '#':
            color = color + (7 - len(color)) * color[len(color) - 1]
        label_to_BC_sort[label] = color

    # apply gradient to labels
    for lbl in ylbls:
        local_label = (lbl.get_text().split(': ('))[0]
        lbl.set_color(label_to_BC_sort[local_label])

    plt.savefig("./result/" + mir_name + "scipy_dendrogram.png", dpi=300)
