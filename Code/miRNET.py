import networkx as nx  # version == 2.3
import pandas as pd
import requests  # HTTP Client for Python
import json  # Standard JSON library
from py2cytoscape import util as cy
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import numpy as np
from colour import Color
import warnings


#   load interactome from String and convert to nx Graph

class MainNet:

    def __init__(self, interactome=None):
        if not isinstance(interactome, pd.DataFrame):
            string = pd.read_csv('./baseData/String_interactome.csv')
            self.G = nx.from_pandas_edgelist(string, 'Source', 'Target')
        else:
            interactome = pd.DataFrame(interactome)
            assert interactome.shape[1] == 2, 'It takes two columns: "Source" and "Target"'
            assert sum(interactome.columns == ['Source', 'Target']) == 2, 'Columns names are not "Source" and "Target"'
            self.G = nx.from_pandas_edgelist(interactome, 'Source', 'Target')
            print('interactome contain ', interactome.shape[0], ' rows')

    def get_LCCnd_centrality(self):
        """
        return: sorted dict of node centrality
        """
        centrality_node = nx.betweenness_centrality(self.get_LCC())
        degree_centrality = nx.degree_centrality(self.get_LCC())
        for k, v in centrality_node.items():
            centrality_node[k] = centrality_node[k] + degree_centrality[k]
        centrality_node = {k: v for k, v in sorted(centrality_node.items(), key=lambda item: item[1], reverse=True)}

        for nods, dct in self.LCC.nodes(data=True):
            dct['BtweenCentrl'] = centrality_node[nods]

        return centrality_node

    def get_LCC(self):
        """
        return: the Largest Connected Component, as NetrworkX object
        """
        CC_G = sorted(nx.connected_component_subgraphs(self.G), key=len, reverse=True)
        self.LCC = CC_G[0]

        # print(len(CC_G), 'total connected components')
        # print(len(CC_G[0].nodes), 'LCC cardinality')

        return CC_G[0]

    def select_nodes(self, gene_set):
        """
        param gene_set: tissue or another gene set
        """

        self.G = self.G.subgraph(gene_set)
        if self.LCC:
            # self.LCC = self.LCC.subgraph(gene_set)
            sub_LCC = sorted(nx.connected_component_subgraphs(self.LCC.subgraph(gene_set)), key=len, reverse=True)
            if len(sub_LCC) > 0:
                self.LCC = sub_LCC[0]
                return 0
            else:
                warnings.warn('!!! len(LCC) == 0 !!! ')
                return 1


def tissue_selector(ans=None, tissue_id=None):
    if ans is None:
        ans = int(str(input('"Human Protein Atlas"(0) or "GTEx"(1) ? ')))
    if ans == 1:
        dt = pd.read_csv(
            './addData/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',
            sep='\t')  # загрузка med(TPM) from GTEx

        print('Gene universe is...')

        labels = sorted(dt.columns)
        index = [i for i in range(0, len(dt.columns))]

        if tissue_id is None:
            for i in index:
                print(index[i], '-----', labels[i])
            print("give me int, or input 'all'", sep='\n')
            tissue_id = str(input('Your choice: '))
        else:
            print(labels[int(tissue_id)], 'was used', sep=' ')

        if tissue_id != 'all':
            tissue = labels[int(tissue_id)]
            tissue_genes = set(dt['Description'][dt[tissue] > 0])
            print('your tissue is ', tissue, ' number of genes: ', len(tissue_genes))

            return tissue_genes
        else:
            return 'all'

    elif ans == 0:
        dt = pd.read_csv(
            './addData/Hum_Atas_normal_tissue.tsv',
            sep='\t')
        print('Gene universe is...')

        labels = sorted(list(map(str, set(dt['Tissue']))))
        index = [i for i in range(0, len(labels))]
        if tissue_id is None:
            for i in index:
                print(index[i], '-----', labels[i])
            print("give me int, or input 'all'", sep='\n')
            tissue_id = str(input('Your choice: '))
        else:
            print(labels[int(tissue_id)], 'was used', sep=' ')

        if tissue_id != 'all':
            tissue = labels[int(tissue_id)]
            tissue_genes = set(dt[(dt['Tissue'] == tissue) & (dt['Level'] != 'Not detected')]['Gene name'])
            print('your tissue is ', tissue, ' number of genes: ', len(tissue_genes))
            return tissue_genes
        else:
            return 'all'
    else:
        return tissue_selector()


class KeyNodesExtractor:

    def __init__(self, MirNet):
        self.key_nodes = dict()
        self.net = MirNet.get_LCC()
        self.node_centrality = MirNet.get_LCCnd_centrality()
        self.graph_features = {'card_LCC': [len(self.net.nodes())],
                               'n_CC': [len(list(nx.connected_component_subgraphs(self.net)))],
                               'transitivity': [nx.transitivity(self.net)],
                               'sh_path': [nx.average_shortest_path_length(self.net) / len(self.net.nodes())]}

    @staticmethod
    def inflection_finder(card_LCC, n_CC, sigma, i=1):
        """
        :param sigma: smoothing
        :param card_LCC: cardinality of the LCC
        :param n_CC: number of connected components in the Network
        :return: the index of the last key node, after the removal of which the network stops rapidly falling apart
        """
        # print('card_LCC:', card_LCC)

        y = gaussian_filter1d(card_LCC, sigma=sigma)
        dy = np.diff(y)  # first derivative

        if dy.size == 0:
            warnings.warn("Can't find any key genes")
            return 'end'

        idx_max_dy = np.argmax(dy)
        # print('ind:', idx_max_dy)
        # print('cardLCC:', card_LCC[idx_max_dy])
        # print('n_CC:', n_CC[idx_max_dy])
        # print('i:', i)
        # print(i)

        if i < 500:     # чтобы пробки рекурсией не выбивало
            if card_LCC[idx_max_dy] > n_CC[idx_max_dy]:
                i += 1
                return KeyNodesExtractor.inflection_finder(card_LCC, n_CC, sigma + 0.2, i)
            else:
                return idx_max_dy

        else:
            ind_dy_zero = [ind for ind, x in enumerate(dy == 0) if x]  # list ind где производная == 0
            # print('ind_dy_zero:', ind_dy_zero)
            ind_dy_zero.append('end')
            for ind in ind_dy_zero:
                if ind == 'end':
                    warnings.warn("Can't find any key genes")
                    return 'end'
                if card_LCC[ind] > n_CC[ind]:
                    continue
                else:
                    return ind

    def extraction(self):

        for k, v in self.node_centrality.items():

            if v == 0:
                break
            self.net.remove_node(k)
            if len(self.net.nodes()) == 0:
                break
            LCC_curent = sorted(nx.connected_component_subgraphs(self.net), key=len, reverse=True)[0]  # get LCC
            # print('LCC_curent:', LCC_curent.nodes())
            self.graph_features['card_LCC'].append(len(LCC_curent.nodes()))
            self.graph_features['n_CC'].append(len(list(nx.connected_component_subgraphs(self.net))))
            self.graph_features['transitivity'].append(nx.transitivity(LCC_curent))
            self.graph_features['sh_path'].append(nx.average_shortest_path_length(LCC_curent) / len(LCC_curent.nodes()))

        # find inflection point of function

        idx_max_dy = KeyNodesExtractor.inflection_finder(card_LCC=self.graph_features['card_LCC'],
                                                         n_CC=self.graph_features['n_CC'],
                                                         sigma=0.0001)
        if idx_max_dy == 'end':
            return 1
        self.graph_features['cutoff_point'] = idx_max_dy

        for i in range(0, idx_max_dy + 1):
            nods = list(self.node_centrality.keys())[i]
            self.key_nodes[nods] = self.node_centrality[nods]

        return self.key_nodes.keys()


class Targets:
    """
    The class extracts and stores the targets of one microRNA and its name
    """

    def __init__(self, path_to_miRTarBase):
        self.miR_dict = {}
        with open(path_to_miRTarBase) as interact:  # import targets from miRTarBase
            for line in interact:
                (key, val) = line.strip().split(';')
                if key in self.miR_dict:
                    self.miR_dict[key] += [val]
                else:
                    self.miR_dict[key] = [val]

    def get_targets(self, miR_name, concat=True):
        """
        miRTarBase содержит разные формы miRNA, так например, miR-21- может соответствовать
        miR-21-3p и miR-21-5p. Функция конкатенирует таргеты всех форм микроРНК и удалаяет дубли,
        которые возникают из-за того, что одной и тойже мишени может соответствовать несколько строк из-за
        разные методов её подтверждения.
        """

        miR_names = self.miR_dict.keys()
        miR_dict = self.miR_dict

        res = set()
        mir_name_app = []
        
        if concat:
            for name in miR_names:
                if miR_name in name:
                    res.update(miR_dict[name])
                    mir_name_app.append(name)
        else:
            if miR_name in miR_dict.keys():
                res = set(miR_dict[miR_name])
                mir_name_app.append(miR_name)

        if not mir_name_app:
            print('miRNA', '"{}"'.format(miR_name), 'not found, use another name')
            return 1

        print('I found a miRNA with name:', *mir_name_app)
        print('and ', len(res), 'unique targets')

        return res

#   visualisation tools

class Plots:

    def __init__(self, MirNet, KeyNodesExtractor, miR_name):
        self.miR_G = MirNet.LCC
        self.key_nodes = KeyNodesExtractor.key_nodes
        self.miR_name = miR_name
        self.centrality_node = MirNet.get_LCCnd_centrality()
        self.card_LCC = KeyNodesExtractor.graph_features['card_LCC']
        self.n_CC = KeyNodesExtractor.graph_features['n_CC']
        self.idx_max_dy = KeyNodesExtractor.graph_features['cutoff_point']

    def central_distr(self):
        miR_G = self.miR_G
        key_nodes = self.key_nodes
        mir_name = self.miR_name

        """

        :param miR_G: a Graph of miRNA
        :param key_nodes: key_nodes from MainMiRNetwork.key_nodes
        :param mir_name: miR name for saving file
        :return: visualisation hist of centrality distribution
        """

        fig = plt.figure()
        ax = fig.add_subplot()
        N, bins, patches = ax.hist([miR_G.nodes[node]['BtweenCentrl'] for node in miR_G.nodes])
        ax.set(xlabel='Centrality',
               ylabel='Count of nodes')
        ax.yaxis.label.set_size(30)
        ax.xaxis.label.set_size(30)
        ax.tick_params(labelsize=20)

        ax.axvline(key_nodes[list(key_nodes.keys())[-1]], color='red', lw='4')

        right_side = ax.spines["right"]
        right_side.set_visible(False)
        top_side = ax.spines["top"]
        top_side.set_visible(False)

        # create gradient (grey_to_red hist path
        grey = Color('#cccccc')
        colors = list(grey.range_to(Color("red"), len(bins) - 1))
        for tmp_color, tmp_patch in zip(colors, patches):
            color = str(tmp_color)
            if len(color) < 7 and color[0] == '#':
                color = color + (7 - len(color)) * color[len(color) - 1]
            tmp_patch.set_facecolor(color)

        fig.set_figwidth(8.5)
        fig.set_figheight(8.5)

        plt.tight_layout()

        plt.savefig('./result/' + mir_name + '_centrality_distr.png', dpi=250)

    def graph_to_cytoscape(self):

        miR_G = nx.Graph(self.miR_G)  # unfreezing of the graph
        centrality_node = self.centrality_node
        rem_CC = 0

        for CC in list(nx.connected_components(miR_G)):
            if len(CC) < 3:
                miR_G.remove_nodes_from(CC)
                rem_CC += 1

        print(rem_CC, ' network components with less than two nodes have been removed', end='\n')

        PORT_NUMBER = 1234
        IP = 'localhost'
        BASE = 'http://' + IP + ':' + str(PORT_NUMBER) + '/v1/'

        requests.delete(BASE + 'session')  # Delete all networks in current session

        cytoscape_network = cy.from_networkx(miR_G)
        cytoscape_network['data']['name'] = 'miR_Net'
        res1 = requests.post(BASE + 'networks', data=json.dumps(cytoscape_network))
        res1_dict = res1.json()
        new_suid = res1_dict['networkSUID']
        requests.get(BASE + 'apply/layouts/force-directed/' + str(new_suid))

        # load and apply style

        res = requests.get(BASE + 'styles/miR_Net_Styles')
        if res.status_code != 200:

            with open('./options/cytoscape_styles/miR_Net_Styles.json') as json_file:
                miR_Net_Styles = json.load(json_file)

            for mapings in range(0, len(miR_Net_Styles['mappings'])):
                if miR_Net_Styles['mappings'][mapings]['visualProperty'] == 'NODE_LABEL_FONT_SIZE':
                    miR_Net_Styles['mappings'][mapings]['points'][1]['value'] = max(centrality_node.values())
                if miR_Net_Styles['mappings'][mapings]['visualProperty'] == 'NODE_SIZE':
                    miR_Net_Styles['mappings'][mapings]['points'][1]['value'] = max(centrality_node.values())
                if miR_Net_Styles['mappings'][mapings]['visualProperty'] == 'NODE_FILL_COLOR':
                    miR_Net_Styles['mappings'][mapings]['points'][1]['value'] = max(centrality_node.values())

            # Create new Visual Style
            res = requests.post(BASE + "styles", data=json.dumps(miR_Net_Styles))

        # Apply it to current network

        requests.get(
            BASE + 'apply/styles/' + 'miR_Net_Styles' + '/' + str(new_suid))  # !Это говно почему-то не работает

    def key_nodes_extractor(self):

        card_LCC = self.card_LCC
        n_CC = self.n_CC
        idx_max_dy = self.idx_max_dy
        mir_name = self.miR_name
        """

        :param card_LCC:
        :param n_CC:
        :param idx_max_dy:
        :param mir_name:
        :return: visualisation plot of key nodes selection
        """

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(card_LCC, linewidth=4, label='Cardinality of the LCC')
        ax.plot(n_CC, linewidth=4, label='Count of CC', color='tab:green', linestyle='dashed')
        # ax.plot(idx_max_dy, card_LCC[idx_max_dy], marker='o', markersize=20, color="red")
        ax.axvline(idx_max_dy, color='red', lw='4')
        ax.minorticks_on()
        ax.grid(which='major',
                color='w',
                linewidth=1.3)
        ax.grid(which='minor',
                color='w',
                linestyle=':')
        ax.set(xlabel='Number of top nodes removed',
               ylabel='LCC cardinality / Count of CC')
        ax.legend()

        right_side = ax.spines["right"]
        right_side.set_visible(False)
        top_side = ax.spines["top"]
        top_side.set_visible(False)

        #    plt.rc('font', size=10)  # controls default text sizes
        plt.rc('axes', labelsize=30)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=20)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=20)  # fontsize of the tick labels
        plt.rc('legend', fontsize=20)  # legend fontsize
        fig.set_figwidth(8)
        fig.set_figheight(8)
        plt.tight_layout()

        plt.savefig('./result/' + mir_name + 'key_nodes_selection.png', dpi=300)
