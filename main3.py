import math
import random
import statistics

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


def centrality_iteration(network, n, centrality_method, iterations):
    dictionary = create_empty_node_dict(network)
    for i in range(iterations):
        print(i)
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
    print('done')
    return dictionary


def avg_std(diction):
    average_dict = {}
    std_dev = {}
    for k, v in diction.items():
        average_dict[k] = sum(v) / len(v)
        std_dev[k] = statistics.pstdev(v)
    return average_dict, std_dev


def t_test(wildtype, disease, std, iterations):
    t_values = {}
    for gene, avg in disease.items():
        if(std[gene]) == 0:
            print(gene)
            t_values[gene] = 0.05
        else:
            t_values[gene] = (avg - wildtype[gene]) / (std[gene] / math.sqrt(iterations))
    return t_values


def avg_difference(wildtype, disease):
    avg_diff = {}
    for gene, avg in disease.items():
        gene_sum = 0
        for metric in avg:
            gene_sum = gene_sum + abs(wildtype[gene] - metric)
        avg_diff[gene] = gene_sum / len(avg)
    return avg_diff


def is_diff_larger(average, wildtype, disease):
    diff_tracker = {}
    true_count = 0
    for gene, metric in disease.items():
        if abs(wildtype[gene] - metric) > average[gene]:
            diff_tracker[gene] = True
            true_count += 1
        else:
            diff_tracker[gene] = False
    print(true_count)
    return diff_tracker


def difference(wildtype, disease):
    diff_tracker = {}
    for gene, metric in disease.items():
        diff_tracker[gene] = abs(metric - wildtype[gene])
    return diff_tracker


def convert_dicts_to_2arr(wt, dis):
    x = []
    y = []
    for gene, metric in dis.items():
        x.append(metric)
        y.append(wt[gene])
    return x, y


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
# # cn_centrality = [nx.degree_centrality(cancer_network), nx.closeness_centrality(cancer_network), nx.betweenness_centrality(cancer_network), nx.pagerank(cancer_network)]
#
# # print(fn_centrality)
# # print(cn_centrality)
#
# fn_centrality_df = pd.DataFrame(fn_centrality).transpose().rename(columns={0: 'Degree', 1: 'Closeness', 2: 'Betweenness', 3: 'PageRank'}).sort_index()
# pretty_print(fn_centrality_df)

# cn_centrality_df = pd.DataFrame(cn_centrality).transpose().rename(columns={0: 'Degree', 1: 'Closeness', 2: 'Betweenness', 3: 'PageRank'})
# print(cn_centrality_df)
#
# print(centrality_iteration(full_network, 17, 'degree'))
dictionary = centrality_iteration(full_network, 17, 'pagerank', 100)
print(dictionary)
# print(avg_std(dictionary)[1])
#
# std_avg = avg_std(dictionary)
# # print(min(list(std_avg.values())))
#
# print(t_test(nx.degree_centrality(full_network), std_avg[0], std_avg[1], 1000))
# print(nx.group_betweenness_centrality(full_network, cancer_genes_in_network))

# cancerous_centrality = fn_centrality_df.iloc[[116252]]
# print(cancerous_centrality)
wt_central = nx.pagerank(full_network)

# print(avg_difference(wt_central, dictionary))
avg_diff_wt_dis = avg_difference(wt_central, dictionary)
print(avg_diff_wt_dis)

cn_central = nx.pagerank(cancer_network)
print(cn_central)
cn_central_list = list(cn_central.values())
print(wt_central)
wt_central_list = list(wt_central.values())
print(is_diff_larger(avg_diff_wt_dis, wt_central, cn_central))

print(len(cn_central))
# plt.hist(avg_diff_wt_dis)
# plt.show()
#
# plt.hist(avg_difference(wt_central, cn_central))
# plt.show()
x, y = convert_dicts_to_2arr(avg_diff_wt_dis, difference(wt_central, cn_central))
fig, ax = plt.subplots()

ax.hist(np.subtract(x, y))
ax.plot([0, 0], [0, len(y)])
# ax.scatter(x, y)
#
# lims = [
#     np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
#     np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
# ]
#
# # now plot both limits against eachother
# ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
# ax.set_aspect('equal')
# ax.set_xlim(lims)
# ax.set_ylim(lims)
# ax.set_xlabel('Difference in wt and Glioblastoma')
# ax.set_ylabel('Average Difference in wt and Random Sampled')
# ax.set_title('Differences in Empirical and Random Sampling Centrality (Degree)')
ax.set_xlabel('Difference between wt and Random Sampling and wt and Glioblastoma')
ax.set_ylabel('Number of Nodes')
ax.set_title('Differences in Empirical and Random Sampling Centrality (PageRank)')
# plt.scatter(avg_diff_wt_dis, difference(wt_central, cn_central))
# print(avg_diff_wt_dis)
# print(difference(wt_central, cn_central))
plt.show()

