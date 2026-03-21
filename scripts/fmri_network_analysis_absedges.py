#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fmri_network_analysis_absedges.py

Purpose:
    Computes graph-theory network metrics from subject-level resting-state
    functional connectivity (FC) matrices using ABSOLUTE edge weights
    (i.e., abs(correlation matrix), diagonal set to zero).
    This is the dissertation version of the network analysis pipeline.

Inputs:
    - Subject-level average FC matrices (.txt), one per subject
    - Cole-Anticevic (CA) network partition label key (.txt):
        <file_keyinfo>
    - Subject list CSV with column 'src_subject_id' 

Outputs (written per subject to <output_dir>/<src_subject_id>/):
    - <src_subject_id>_nodal_stats.csv   : nodal-level graph metrics
    - <src_subject_id>_distance.csv      : pairwise edge distances

Output (written to <output_dir>/):
    - network_summary_stats_<first_subj>_to_<last_subj>.csv : group summary

Metrics computed:
    Nodal:   strength, normalized strength, closeness centrality,
             betweenness centrality, eigenvector centrality, PageRank,
             clustering coefficient, participation coefficient (CA + Louvain)
    Global:  mean degree, average shortest path length,
             modularity (CA partition + Louvain), number of Louvain communities,
             average clustering coefficient, mean participation coefficient

Dependencies:
    numpy, pandas, networkx, community (python-louvain)

Usage:
    Update the USER CONFIGURATION section below, then run:
        python fmri_network_analysis_absedges.py
    Or submit via SLURM:
        sbatch run_fmri_network_analysis.slurm
"""

# Note: on SBU Seawulf cluster, the pandas hashtable module is in a
# non-standard location and must be imported explicitly before pandas loads.
# Remove or comment out this line if running on other systems.
# from pandas._libs import hashtable

import numpy as np
import pandas as pd
import os
import networkx as nx
import community  # python-louvain


# ==============================================================================
# USER CONFIGURATION — update these paths before running
# ==============================================================================

# Row range of subjects to process from the subject list (0-based, exclusive end)
# Adjust to distribute subjects across jobs
row_range = np.arange(0, 725)

# Cole-Anticevic network partition label key file
file_keyinfo = "/path/to/script/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"

# Directory containing subject average FC matrix files
input_dir = "/path/to/HCA_stat/subj_Avecormat_txt"

# Output directory (will contain one subfolder per subject)
output_dir = "/path/to/HCA_stat/network_output"

# Subject list CSV
subjfile = "/path/to/script/HCA_subject_list.csv"

# ==============================================================================


# ------------------------------------------------------------------------------
# Load partition metadata
# ------------------------------------------------------------------------------

keyinfo = pd.read_csv(file_keyinfo, sep="\t")

# Parcel labels and CA network partition labels
parcelLabels    = keyinfo['LABEL']
partitionLabels = keyinfo['NETWORK']


# ------------------------------------------------------------------------------
# Load subject list
# ------------------------------------------------------------------------------

subjinfo = pd.read_csv(subjfile)
subj_list = subjinfo.loc[row_range, 'src_subject_id']


# ------------------------------------------------------------------------------
# Helper function: participation coefficient (weighted graph)
# ------------------------------------------------------------------------------

def participation_coefficient(G, module_partition):
    """
    Compute the participation coefficient of each node in G.

    Based on Guimera & Amaral (2005). Adapted for weighted graphs:
    within-module degree is the sum of weighted edges to nodes in the
    same module, and total degree is the weighted node strength.

    Parameters
    ----------
    G : networkx.Graph
        Weighted undirected graph.
    module_partition : dict
        Maps each module label to a list of node indices in that module.

    Returns
    -------
    dict
        Maps each node index to its participation coefficient.
    """
    pc_dict = {}

    for m in module_partition.keys():
        M = set(module_partition[m])
        for v in M:
            # Total weighted degree (strength) of node v
            degree = float(nx.degree(G=G, nbunch=v, weight='weight'))
            # Within-module weighted degree: sum of edge weights to nodes in M
            wm_degree = float(sum([G.edges[u, v]['weight'] for u in M if (u, v) in G.edges()]))
            # Participation coefficient: 1 - (within-module fraction)^2
            pc_dict[v] = 1 - ((wm_degree / degree) ** 2)

    return pc_dict


# ------------------------------------------------------------------------------
# Initialize group-level summary DataFrame
# ------------------------------------------------------------------------------

df_sumStat = pd.DataFrame(columns=[
    'src_subject_id',
    'mean_degree',
    'average_shortest_path_length',
    'modularity_CAnet',
    'modularity_Louvain',
    'Louvain_community',
    'average_clustering_coef',
    'mean_participation_coef_CAnet',
    'mean_participation_coef_Louvain'
])


# ==============================================================================
# MAIN LOOP — process each subject
# ==============================================================================

for subject in subj_list:
    subj_src_id = subject
    # Extract numeric ID: FC matrices are named with numeric ID only
    subj_num_id = subj_src_id.partition('HCA')[2].split('_')[0]
    subj_out_dir = os.path.join(output_dir, subj_src_id)

    # ── Load FC matrix ────────────────────────────────────────────────────────

    matrix_f   = subj_num_id + '_rest_CormatAveRuns.txt'
    subj_matrix = os.path.join(input_dir, matrix_f)

    # Skip subject if FC matrix file is not found
    if not os.path.exists(subj_matrix):
        continue

    # Create subject output directory if needed
    if not os.path.exists(subj_out_dir):
        os.mkdir(subj_out_dir)

    # ── Prepare matrix: set diagonal to 0, absolutise ────────────────────────

    matrix = np.genfromtxt(subj_matrix)

    # The Fisher r-to-z transform (wb_command) outputs diagonal values ~7.25;
    # set diagonal to 0 before constructing the graph
    matrixdiag0 = matrix.copy()
    np.fill_diagonal(matrixdiag0, 0)

    # Absolutise all edge weights
    matrixdiag0 = abs(matrixdiag0)

    # ── Build graph ───────────────────────────────────────────────────────────

    G = nx.from_numpy_array(matrixdiag0)
    G.remove_edges_from(list(nx.selfloop_edges(G)))

    # ── Nodal degree (weighted strength) ─────────────────────────────────────

    strength = G.degree(weight='weight')
    strengths = {node: val for (node, val) in strength}
    nx.set_node_attributes(G, strengths, 'strength')

    # Normalized strength: divide by (N - 1)
    normstrengths = {node: val * 1 / (len(G.nodes) - 1) for (node, val) in strength}
    nx.set_node_attributes(G, normstrengths, 'strengthnorm')

    normstrengthlist = np.array([val * 1 / (len(G.nodes) - 1) for (node, val) in strength])
    mean_degree = np.sum(normstrengthlist) / len(G.nodes)

    # ── Distance attribute (inverse of weight, for path-based metrics) ────────

    # For correlation-based networks, distance = 1 / |weight|
    # so that strongly correlated nodes are "closer"
    G_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in G.edges(data='weight')}
    nx.set_edge_attributes(G, G_distance_dict, 'distance')

    # ── Centrality metrics ────────────────────────────────────────────────────

    closeness   = nx.closeness_centrality(G, distance='distance')
    nx.set_node_attributes(G, closeness, 'closecent')

    betweenness = nx.betweenness_centrality(G, weight='distance', normalized=True)
    nx.set_node_attributes(G, betweenness, 'betweencent')

    eigen       = nx.eigenvector_centrality(G, weight='weight')
    nx.set_node_attributes(G, eigen, 'eigen')

    pagerank    = nx.pagerank(G, weight='weight')
    nx.set_node_attributes(G, pagerank, 'pg')

    # ── Average shortest path length ──────────────────────────────────────────

    ave_shortest_path = nx.average_shortest_path_length(G, weight='distance')

    # ── Modularity ────────────────────────────────────────────────────────────

    # Modularity under the CA prior partition
    # Node index must start from 0; CA label key is 1-indexed so subtract 1
    prior_part = dict(zip(keyinfo['INDEX'] - 1, partitionLabels))
    prior_part_modularity = community.modularity(prior_part, G, weight='weight')

    # Louvain community detection and modularity
    part = community.best_partition(G, weight='weight')
    part_modularity = community.modularity(part, G, weight='weight')

    # ── Participation coefficient ─────────────────────────────────────────────

    # Invert partition dict: module label → list of node indices
    ca_module_partition = {}
    for n, m in prior_part.items():
        ca_module_partition.setdefault(m, []).append(n)

    louvain_module_partition = {}
    for n, m in part.items():
        louvain_module_partition.setdefault(m, []).append(n)

    ca_part_participationcoef      = participation_coefficient(G, ca_module_partition)
    louvain_part_participationcoef = participation_coefficient(G, louvain_module_partition)

    # Convert to sorted DataFrames (sorted by node index to align with nodal df)
    df_ca_parcoef = pd.DataFrame.from_dict(ca_part_participationcoef, orient='index').sort_index()
    df_ca_parcoef.columns = ['participationcoef_CAnet']

    df_louvain_parcoef = pd.DataFrame.from_dict(louvain_part_participationcoef, orient='index').sort_index()
    df_louvain_parcoef.columns = ['participationcoef_Louvain']

    # Mean participation coefficient across all nodes (Chan et al., 2014)
    ca_mean_partcoef      = df_ca_parcoef['participationcoef_CAnet'].mean()
    louvain_mean_partcoef = df_louvain_parcoef['participationcoef_Louvain'].mean()

    # ── Clustering coefficient ────────────────────────────────────────────────

    clustering = nx.clustering(G, weight='weight')
    nx.set_node_attributes(G, clustering, 'cc')

    ave_clustering = nx.average_clustering(G, weight='weight')

    # ── Build and save nodal-level output ─────────────────────────────────────

    df_node = pd.DataFrame({
        'strenth':           strengths.values(),
        'normstrength':      normstrengths.values(),
        'closecent':         closeness.values(),
        'betweencent':       betweenness.values(),
        'eigencent':         eigen.values(),
        'pgcent':            pagerank.values(),
        'NETWORKKEY_LOUVAIN': part.values(),
        'clustercoef':       clustering.values()
    })

    # Shift Louvain community index to start from 1 (consistent with CA partition)
    df_node['NETWORKKEY_LOUVAIN'] = df_node['NETWORKKEY_LOUVAIN'] + 1

    df_node = df_node.join(df_ca_parcoef['participationcoef_CAnet'])
    df_node = df_node.join(df_louvain_parcoef['participationcoef_Louvain'])

    # Prepend parcel INDEX column (1-based, to match CA partition key)
    df_node.insert(0, 'INDEX', np.arange(1, len(df_node) + 1))
    df_node.index = np.arange(1, len(df_node) + 1)

    # Pairwise edge distances
    df_distance = pd.DataFrame(list(G.edges(data='distance')),
                               columns=['Source', 'Target', 'Distance'])
    df_distance.index = np.arange(1, len(df_distance) + 1)

    # Save to CSV
    node_f     = os.path.join(subj_out_dir, subj_src_id + '_nodal_stats.csv')
    distance_f = os.path.join(subj_out_dir, subj_src_id + '_distance.csv')

    df_node.to_csv(node_f,     sep=',', header=True, index=False)
    df_distance.to_csv(distance_f, sep=',', header=True, index=False)

    # ── Append to group summary ───────────────────────────────────────────────

    temp_df = pd.DataFrame([[
        subj_src_id,
        mean_degree,
        ave_shortest_path,
        prior_part_modularity,
        part_modularity,
        len(set(part.values()).union()),
        ave_clustering,
        ca_mean_partcoef,
        louvain_mean_partcoef
    ]], columns=[
        'src_subject_id',
        'mean_degree',
        'average_shortest_path_length',
        'modularity_CAnet',
        'modularity_Louvain',
        'Louvain_community',
        'average_clustering_coef',
        'mean_participation_coef_CAnet',
        'mean_participation_coef_Louvain'
    ])

    df_sumStat = pd.concat([df_sumStat, temp_df], axis=0, ignore_index=True)


# ==============================================================================
# Save group-level summary
# ==============================================================================

df_sumStat.index = np.arange(1, len(df_sumStat) + 1)

sumStat_filename = (
    'network_summary_stats_'
    + subj_list.iloc[0]
    + '_to_'
    + subj_list.iloc[-1]
    + '.csv'
)

df_sumStat.to_csv(os.path.join(output_dir, sumStat_filename),
                  sep=',', header=True, index=False)
