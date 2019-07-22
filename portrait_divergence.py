#!/usr/bin/env python
# -*- coding: utf-8 -*-

# portrait_divergence.py
# Jim Bagrow
# Last Modified: 2018-04-24

import sys, os
import tempfile
import argparse
from collections import Counter
import numpy as np
import networkx as nx
from scipy.stats import entropy


def portrait_cpp(graph, fname=None, keepfile=False):
    """Compute and generate portrait of graph using compiled B_matrix
    executable.
    
    Return matrix B where B[i,j] is the number of starting nodes in graph with
    j nodes in shell i
    """
    # file to save to:
    f = fname
    if fname is None:
        f = next(tempfile._get_candidate_names())
    
    # make sure nodes are 0,...,N-1 integers:
    graph = nx.convert_node_labels_to_integers(graph)
    
    # write edgelist:
    nx.write_edgelist(graph, f+".edgelist", data=False)
    
    # make B-matrix:
    os.system("./B_matrix {}.edgelist {}.Bmat > /dev/null".format(f, f))
    portrait = np.loadtxt("{}.Bmat".format(f))
    
    # clean up:
    if not keepfile:
        os.remove(f+".edgelist")
        os.remove(f+".Bmat")
    
    return portrait


def portrait_py(graph):
    """Return matrix B where B[i,j] is the number of starting nodes in graph
    with j nodes in shell i.
    
    If this function is too slow, consider portrait_cpp() instead.
    """
    dia = 500 #nx.diameter(graph)
    N = graph.number_of_nodes()
    # B indices are 0...dia x 0...N-1:
    B = np.zeros((dia+1,N)) 
    
    max_path = 1
    adj = graph.adj
    for starting_node in graph.nodes():
        nodes_visited = {starting_node:0}
        search_queue = [starting_node]
        d = 1
        while search_queue:
            next_depth = []
            extend = next_depth.extend
            for n in search_queue:
                l = [i for i in adj[n] if i not in nodes_visited] 
                extend(l)
                for j in l:
                    nodes_visited[j] = d
            search_queue = next_depth
            d += 1
            
        node_distances = nodes_visited.values()
        max_node_distances = max(node_distances)
        
        curr_max_path = max_node_distances
        if curr_max_path > max_path:
            max_path = curr_max_path
        
        # build individual distribution:
        dict_distribution = dict.fromkeys(node_distances, 0)
        for d in node_distances:
            dict_distribution[d] += 1
        # add individual distribution to matrix:
        for shell,count in dict_distribution.items():
            B[shell][count] += 1
        
        # HACK: count starting nodes that have zero nodes in farther shells
        max_shell = dia
        while max_shell > max_node_distances:
            B[max_shell][0] += 1
            max_shell -= 1
    
    return B[:max_path+1,:]


portrait = portrait_py
#portrait = portrait_cpp


def weighted_portrait(G, paths=None, binedges=None):
    """Compute weighted portrait of G, using Dijkstra's algorithm for finding
    shortest paths. G is a networkx object.
    
    Return matrix B where B[i,j] is the number of starting nodes in graph with
    j nodes at distance d_i <  d < d_{i+1}.
    """
    # all pairs path lengths
    if paths is None:
        paths = list(nx.all_pairs_dijkstra_path_length(G))
    
    if binedges is None:
        unique_path_lengths  = _get_unique_path_lengths(G, paths=paths)
        sampled_path_lengths = np.percentile(unique_path_lengths, np.arange(0, 101, 1))
    else:
        sampled_path_lengths = binedges
    UPL = np.array(sampled_path_lengths)
    
    l_s_v = []
    for i,(s,dist_dict) in enumerate(paths):
        distances = np.array(list(dist_dict.values()))
        s_v,e = np.histogram(distances, bins=UPL)
        l_s_v.append(s_v)
    M = np.array(l_s_v)
    
    B = np.zeros((len(UPL)-1, G.number_of_nodes()+1))
    for i in range(len(UPL)-1):
        col = M[:,i] # ith col = numbers of nodes at d_i <= distance < d_i+1
        for n,c in Counter(col).items():
            B[i,n] += c
    
    return B


def _get_unique_path_lengths(graph, paths=None):
    if paths is None:
        paths = list(nx.all_pairs_dijkstra_path_length(graph))

    unique_path_lengths = set()
    for starting_node,dist_dict in paths:
        unique_path_lengths |= set(dist_dict.values())
    unique_path_lengths = sorted(list(unique_path_lengths))
    return unique_path_lengths


def pad_portraits_to_same_size(B1,B2):
    """Make sure that two matrices are padded with zeros and/or trimmed of
    zeros to be the same dimensions.
    """
    ns,ms = B1.shape
    nl,ml = B2.shape
    
    # Bmats have N columns, find last *occupied* column and trim both down:
    lastcol1 = max(np.nonzero(B1)[1])
    lastcol2 = max(np.nonzero(B2)[1])
    lastcol = max(lastcol1,lastcol2)
    B1 = B1[:,:lastcol+1]
    B2 = B2[:,:lastcol+1]
    
    BigB1 = np.zeros((max(ns,nl), lastcol+1))
    BigB2 = np.zeros((max(ns,nl), lastcol+1))
    
    BigB1[:B1.shape[0],:B1.shape[1]] = B1
    BigB2[:B2.shape[0],:B2.shape[1]] = B2
    
    return BigB1, BigB2


def _graph_or_portrait(X):
    """Check if X is a nx (di)graph. If it is, get its portrait. Otherwise
    assume it's a portrait and just return it.
    """
    if isinstance(X, (nx.Graph, nx.DiGraph)):
        return portrait(X)
    return X


def portrait_divergence(G, H):
    """Compute the network portrait divergence between graphs G and H."""
    
    BG = _graph_or_portrait(G)
    BH = _graph_or_portrait(H)
    BG, BH = pad_portraits_to_same_size(BG,BH)
    
    L, K = BG.shape
    V = np.tile(np.arange(K),(L,1))
    
    XG = BG*V / (BG*V).sum()
    XH = BH*V / (BH*V).sum()
    
    # flatten distribution matrices as arrays:
    P = XG.ravel()
    Q = XH.ravel()
    
    # lastly, get JSD:
    M = 0.5*(P+Q)
    KLDpm = entropy(P, M, base=2)
    KLDqm = entropy(Q, M, base=2)
    JSDpq = 0.5*(KLDpm + KLDqm)
    
    return JSDpq


def portrait_divergence_weighted(G,H, bins=None, binedges=None):
    """Network portrait divergence between two weighted graphs.
    
    bins = width of bins in percentiles
    binedges = vector of bin edges
    bins and binedges are mutually exclusive
    """
    
    # get joint binning:
    paths_G = list(nx.all_pairs_dijkstra_path_length(G))
    paths_H = list(nx.all_pairs_dijkstra_path_length(H))
    
    # get bin_edges in common for G and H:
    if binedges is None:
        if bins is None:
            bins = 1
        UPL_G = set(_get_unique_path_lengths(G, paths=paths_G))
        UPL_H = set(_get_unique_path_lengths(H, paths=paths_H))
        unique_path_lengths = sorted(list(UPL_G | UPL_H))
        binedges = np.percentile(unique_path_lengths, np.arange(0, 101, bins))
    
    # get weighted portraits:
    BG = weighted_portrait(G, paths=paths_G, binedges=binedges)
    BH = weighted_portrait(H, paths=paths_H, binedges=binedges)
    
    return portrait_divergence(BG, BH)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


if __name__ == '__main__':

    description_str="""
Compute network portrait divergences between a pair of networks. Prints the
portrait divergence to STDOUT.
"""
    
    epilog_str = """
Network file format
===================
Files may be edgelists (see networkx.read_edgelist()) or graphml files
(see networkx.read_graphml()).

Edgelists specify the network in a two- or three-column matrix format:
    node_i <delimiter> node_j <delimiter> [weight_ij]
    ...
Columns are delimited with whitespace by default.

""".format(fn=sys.argv[0])
    
    # process arguments
    parser = argparse.ArgumentParser(description=description_str,
            epilog=epilog_str, formatter_class=CustomFormatter)
    parser.add_argument('filename1', metavar='graph1', help='Filename for the first network')
    parser.add_argument('filename2', metavar='graph2', help='Filename for the second network')
    parser.add_argument('-d', '--directed', action='store_true',
            help='treat networks as directed.')
    parser.add_argument('-w', '--weighted', default=False, nargs="?", const='weight',
            help='treat networks as weighted. Optional argument WEIGHTED gives the edge attribute key for the weights')
    parser.add_argument('-b', '--binning', default=1, nargs="?", 
            help='width of portrait bins in percentiles if networks are considered weighted')
    parser.add_argument("-e", "--delimiter", default=None,
            help="specify the column delimiter used for edgelist files. Default: contiguous whitespace. Ignored if --graphml used")
    parser.add_argument('--graphml', action='store_true', 
            help='input files are .graphml file instead of two/three-column edgelists')
    parser.add_argument('--cpp', action='store_true', 
            help='use faster C++ implementation of the network portraits. Requires B_matrix executable to be installed. Does not support directed or weighted graphs.')
    args = parser.parse_args()
    #print(args)
    
    if args.cpp:
        if args.directed or args.weighted:
            sys.exit("The C++ code does not currently support directed or weighted graphs. Use the Python code instead. Exiting...")
        portrait = portrait_cpp
    
    # read networks
    for f in [args.filename1, args.filename2]:
        if not os.path.exists(f):
            sys.exit("File {} does not exist. Exiting...".format(f))
    if args.graphml:
        G = nx.read_graphml(args.filename1)
        H = nx.read_graphml(args.filename2)
    else:
        create_using = nx.Graph
        if args.directed:
            create_using = nx.DiGraph
        data = False
        if args.weighted:
            data = [(args.weighted, float)]
        
        G = nx.read_edgelist(args.filename1, create_using=create_using(), data=data)
        H = nx.read_edgelist(args.filename2, create_using=create_using(), data=data)
        
    if args.weighted is False:
        Djs = portrait_divergence(G, H)
    else:
        Djs = portrait_divergence_weighted(G, H, bins=args.binning)
    
    #print(nx.info(G))
    #print(nx.info(H))
    print(Djs)







