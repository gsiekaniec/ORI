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
    cutoff = 0.01
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

if __name__ == '__main__':
    main()







#N = {}
#D = {}
#for i in range(l):
#    N[i] = [n for n in G.neighbors(i)]
#    D[i] = len(N[i])
#
#def get_N_D(n, G):
#    N = [x for x in G.neighbors(n)]
#    return N, len(N)
#
#coeff = {}
#a = []
#for i in range(l):
#    neigh = N[i]
#    di = D[i]
#    if di in (0,1): continue
#    dd = di * (di-1)
#    tmp_graph = nx.Graph()
#    for n1, n2 in itertools.combinations(neigh, 2):
#        tmp_graph.add_nodes_from([n1, n2])
#        if adj_mat[n1][n2] == 1:
#            tmp_graph.add_edge(n1, n2)
#    nv = tmp_graph.number_of_edges()
#    coeff[i] = nv/dd
#    a.append(coeff[i])
#
#def get_graph(L, G):
#    nG = nx.Graph()
#    for n1, n2 in itertools.combinations(L, 2):
#        nG.add_nodes_from([n1, n2])
#        if adj_mat[n1][n2] == 1:
#            nG.add_edge(n1, n2)
#    return nG
#
#def get_coeff(node, graph):
#    N, D = get_N_D(node, graph)
#    nG = get_graph(N, graph)
#    nb_e = nG.number_of_edges()
#    return 2 * (nb_e / (D * (D-1)))
#
#
#def get_overlap(G, c, all_c):
#    overlap = []
#    for node in c:
#        for cl in all_c:
#            if node in cl:
#                overlap.append(cl)
#    return overlap
#
#mc = list(nx.find_cliques_recursive(G))
#mc = sorted(mc, key=len, reverse=True)
#with open('max_cliques.txt', 'w') as f_out:
#    for cl in mc:
#        if len(cl) > 1:
#            f_out.write(' '.join(list(map(str, cl))) + '\n')
#
#refined_clique = []
#nodes = set() 
#for c in mc:
#    a = 0
#    for i in c:
#        if i in nodes:
#            a = 1
#            break
#    if a == 1: continue
#    refined_clique.append(c)
#    for i in c: nodes.add(i)
#
#print(refined_clique)
#with open('max_cliquesR.txt', 'w') as f_out:
#    for cl in refined_clique:
#        if len(cl) > 1:
#            f_out.write(' '.join(list(map(str, cl))) + '\n')
#    
#
#non_over_cliques = []
#for i in range(len(mc)):
#    c = mc[i]
#    tmp = mc.copy()
#    del tmp[i]
#    o = 0
#    for cls in mc:
#        if any(i in cls for i in c):
#            o += 1
#    if not o: non_over_cliques.append(c)
#print(non_over_cliques)
#
#sys.exit()
#CP = nx.Graph(G)
#for c in mc:
#    overlap = get_overlap(CP, c, mc)
#    #for o in overlap:
#        
#
#
#
#
#
#clique_by_node = {}
#for i in range(l):
#    if i not in clique_by_node:
#        clique_by_node[i] = []
#    for c in mc:
#        if i in c:
#            clique_by_node[i].append(c)
#
#overlap_cliques = []
#for s1, s2 in itertools.combinations(mc, 2):
#    if bool(set(s1) & set(s2)):
#        overlap_cliques.append((s1, s2))
#
#for c1, c2 in overlap_cliques:
#    common = set(c1) & set(c2)
#    g = get_graph(list(common), None)
#    for n1, n2 in itertools.combinations(common, 2):
#        N1, D1 = get_N_D(n1, g)
#        N2, D2 = get_N_D(n2, g)
#        a = len(set(N1).union(set(N2)))
#        den = a * (a-1)
#        l = len(N1)+len(N2)
#        print(str((2*l)/den))


    


#print(overlap_cliques)
#print(list(nx.find_cliques_recursive(G)))