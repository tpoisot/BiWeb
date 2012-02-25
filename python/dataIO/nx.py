from ..bipartite_class import *
import networkx as nx


def nxImport(nxg,name=''):
    if not nx.algorithms.bipartite.is_bipartite(nxg):
        print 'This graph is not bipartite'
        return 0
    else:
        X,Y = nx.algorithms.bipartite.sets(nxg)
        top_nodes = list(X)
        bottom_nodes = list(Y)
        W = np.zeros((len(top_nodes),len(bottom_nodes)))
        einf =[(u,v,d) for (u,v,d) in nxg.edges(data=True)]
        for edge in einf:
            sp1 = edge[0]
            sp2 = edge[1]
            edata = edge[2]
            if 'weight' in edata:
                ls = edata['weight']
            else:
                ls = 1
            if sp1 in top_nodes:
                W[top_nodes.index(sp1)][bottom_nodes.index(sp2)] = ls
            else:
                W[top_nodes.index(sp2)][bottom_nodes.index(sp1)] = ls
        Bip = bipartite(W)
        Bip.name = name
        Bip.upnames = top_nodes
        Bip.lonames = bottom_nodes
        return Bip
