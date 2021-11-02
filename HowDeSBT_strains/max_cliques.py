#!/usr/bin/env python3
# coding: utf-8

import random
import sys
import networkx as nx
import numpy
import itertools
try:
    import pygraphviz
    from networkx.drawing.nx_agraph import write_dot
    print("using package pygraphviz")
except ImportError:
    try:
        import pydot
        from networkx.drawing.nx_pydot import write_dot
        print("using package pydot")
    except ImportError:
        print()
        print("Both pygraphviz and pydot were not found ")
        print("see  https://networkx.github.io/documentation/latest/reference/drawing.html")
        print()
        raise

def get_mcliques(G, std=True):
    if std:
        return sorted(list(nx.find_cliques_recursive(G)), key=len, reverse=True)
    else:
        return list(nx.find_cliques_recursive(G))

def get_overlapping_cliques(mc):
    overlaps = []
    for c1, c2 in itertools.combinations(mc, 2):
        if bool(set(c1) & set(c2)):
            overlaps.append((c1, c2))
    return sorted(overlaps)

def get_clco(i, j, Ni, Nj, G):
    uN = set(Ni).union(set(Nj))
    d = len(uN)
    if d <= 1: return 1
    dd = d * (d-1)
    tmp_g = nx.Graph()
    tmp_g.add_nodes_from(uN)
    for n1, n2 in itertools.combinations(uN, 2):
        if G.has_edge(n1, n2):
            tmp_g.add_edge(n1, n2)
    nbe = tmp_g.number_of_edges()
    return 2 * nbe / dd

def main():
    G = nx.Graph()

    in_path = sys.argv[1]

    adj_mat = []

    with open(in_path, 'r') as f_in:
        for line in f_in:
            line = line.strip()
            adj_mat.append(list(map(int, line.split(' '))))

    l = len(adj_mat)
    G.add_nodes_from(list(range(l)))

    for i in range(l):
        for j in range(l):
            if adj_mat[i][j] == 1:
                if i != j:
                    G.add_edge(i, j)

    write_dot(G, "initial_graph.dot")


    if len(sys.argv) > 2:
        cutoff = float(sys.argv[2])
        mc = get_mcliques(G)
        with open('max_cliques.txt', 'w') as f_out:
            for cl in mc:
                if len(cl) > 1:
                    f_out.write(' '.join(sorted(list(map(str, cl)))) + '\n')

        overlaps = get_overlapping_cliques(mc)
        print(overlaps)
        while overlaps: # while overlapping maximal cliques are found
            c1, c2 = set(overlaps[0][0]), set(overlaps[0][1]) # get first overlap
            intersect = c1.intersection(c2)
            c1_bis, c2_bis = c1.copy(), c2.copy()
            for i in intersect:
                c1_bis.remove(i); c2_bis.remove(i)
            c1n, c2n = random.sample(c1_bis, 1)[0], random.sample(c2_bis, 1)[0] # get one node from each clique
            N1, N2 = G.neighbors(c1n), G.neighbors(c2n) # get neihbors of each node in
            clco = get_clco(c1n, c2n, N1, N2, G) # compute clco
            if clco >= cutoff:
                for n1, n2 in itertools.combinations(c1.union(c2), 2): #merge cliques
                    G.add_edge(n1, n2)
            else: # compute minimum cut to split c1, c2
                tmp_graph = nx.Graph()
                tmp_graph.add_nodes_from(c1.union(c2))
                for c in (c1, c2):
                    for n1, n2 in itertools.combinations(c, 2):
                        tmp_graph.add_edge(n1, n2)
                cuts = nx.minimum_edge_cut(tmp_graph)
                G.remove_edges_from(cuts) # split c1, c2 in global graph
            mc = get_mcliques(G) # recompute maximal cliques
            overlaps = get_overlapping_cliques(mc) # recompute overlapping cliques


        write_dot(G, "non_ov_clique_graph.dot")
        with open('max_cliques_wo_overlap.txt', 'w') as f_out:
            for cl in mc:
                if len(cl) > 1:
                    f_out.write(' '.join(sorted(list(map(str, cl)))) + '\n')
    else:
        connected_components = nx.connected_components(G)
        with open('max_cliques_wo_overlap.txt', 'w') as f_out:
            for cc in connected_components:
                if len(cc) > 1:
                    f_out.write(' '.join(sorted(list(map(str, cc)))) + '\n')


if __name__ == '__main__':
    main()


