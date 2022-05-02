import random

from networkx import NetworkXError, nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def pretty_print(matrix):
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)
    print(pd.DataFrame(matrix))


def centrality_tests(network):
    return [nx.degree_centrality(network), nx.closeness_centrality(network), nx.betweenness_centrality(network), nx.pagerank(network)]


def remove_n_nodes(network, n):
    full_nodes = network.nodes
    random_sample = random.sample(full_nodes, n)
    sampled_network = network.copy()
    sampled_network.remove_nodes_from(random_sample)
    return sampled_network


def create_empty_node_dict(network):
    dictionary = dict()
    for node in network.nodes:
        dictionary[node] = []
    return dictionary


def centrality_iteration(network, n, centrality_method):
    dictionary = create_empty_node_dict(network)
    for i in range(3):
        # centrality = centrality_tests(remove_n_nodes(network, n))
        if centrality_method == 'degree':
            centrality = nx.degree_centrality(remove_n_nodes(network, n))
        elif centrality_method == 'closeness':
            centrality = nx.closeness_centrality(remove_n_nodes(network, n))
        elif centrality_method == 'betweenness':
            centrality = nx.betweenness_centrality(remove_n_nodes(network, n))
        # elif centrality_method == 'pagerank':
        else:
            centrality = nx.pagerank(remove_n_nodes(network, n))
        for node in centrality:
            dictionary[node].append(centrality[node])

    return dictionary


def average_dict_val(dictionary):
    average_dict = {}
    for k, v in dictionary.items():
        average_dict[k] = sum(v) / float(len(v))



full_edgelist = pd.read_csv('full_network.txt', sep='\t', lineterminator='\n')
cancer = pd.read_csv('cancer_genes.txt', sep='\t', lineterminator='\n')
full_edgelist = full_edgelist[['BioGRID ID Interactor A', 'BioGRID ID Interactor B']]
cancer = cancer[['#BIOGRID ID']]

# pretty_print(full_edgelist)

cancer_genes_in_network = []

for index, row in cancer.iterrows():
    if row['#BIOGRID ID'] in set(full_edgelist['BioGRID ID Interactor A']) or row['#BIOGRID ID'] in set(
            full_edgelist['BioGRID ID Interactor B']):
        cancer_genes_in_network.append(row['#BIOGRID ID'])

print(cancer_genes_in_network)
# print(len(cancer_genes_in_network))

full_network = nx.from_pandas_edgelist(full_edgelist, source='BioGRID ID Interactor A',
                                       target='BioGRID ID Interactor B')

cancer_network = nx.from_pandas_edgelist(full_edgelist, source='BioGRID ID Interactor A',
                                         target='BioGRID ID Interactor B')

cancer_network.remove_nodes_from(cancer_genes_in_network)

print(full_network.nodes)
print(cancer_network.nodes)

# print(centrality_tests(full_network))
#
# fn_centrality = [nx.degree_centrality(full_network), nx.closeness_centrality(full_network),
#                  nx.betweenness_centrality(full_network), nx.pagerank(full_network)]
# cn_centrality = [nx.degree_centrality(cancer_network), nx.closeness_centrality(cancer_network), nx.betweenness_centrality(cancer_network), nx.pagerank(cancer_network)]

# print(fn_centrality)
# print(cn_centrality)

# fn_centrality_df = pd.DataFrame(fn_centrality).transpose().rename(columns={0: 'Degree', 1: 'Closeness', 2: 'Betweenness', 3: 'PageRank'})
# print(fn_centrality_df)
#
# cn_centrality_df = pd.DataFrame(cn_centrality).transpose().rename(columns={0: 'Degree', 1: 'Closeness', 2: 'Betweenness', 3: 'PageRank'})
# print(cn_centrality_df)
#
# print(centrality_iteration(full_network, 17, 'degree'))
# print(average_dict_val(centrality_iteration(full_network, 17, 'degree')))

print(nx.group_betweenness_centrality(full_network, cancer_genes_in_network))
