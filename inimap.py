'''This module constructs weighted graph and topgraph initial mappings ''' 
import networkx as nx
from FiDLS.utils import graph_of_circuit
from FiDLS.utils import hub
from vfs import Vf
'''Always put local parameters before global ones''' 

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
    
def is_embeddable(g, H, anchor, stop):
    '''check if a small graph g is embeddable in a large H, anchor is bool
        g, H (Graph)
        anchor (bool): whether or not mapping anchor of g to that of H
        stop (float): time limit for vf2
    '''
    vf2 = Vf()
    result = {} 
    if anchor: result[hub(g)] = hub(H)
    result = vf2.dfsMatch(g, H, result, stop)
    lng = len(nx.nodes(g))
    if len(result) == lng:
        return True, result   
    return False, result

'''The edge by edge search in 'y' can be sped up if we consider in a bipartite way'''
def top_z_P(i, L_temp, C, G, anchor, stop):
    '''L_temp is the list of the indices of the first cnot gates corresponding to edges in g_of_c'''
    g = nx.Graph()
    for s in L_temp[:i+1]:
        g.add_edge(C[s][0], C[s][1])
    test = is_embeddable(g, G, anchor, stop)
    return test[0]

def search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop):
    '''L_temp is the list of the indices of the first cnot gates corresponding to edges in g_of_c'''
    if not type(test_number) == int: raise Exception('Only consider integers')
    if no_bound == yes_bound + 1: return yes_bound
    if top_z_P(test_number, L_temp, C, G, anchor, stop):
        yes_bound = test_number
        if test_number == no_bound: return test_number
        test_number = yes_bound + max(1, (no_bound - yes_bound)//2)
        return search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop)
    else:
        if test_number == yes_bound + 1: return yes_bound
        no_bound = test_number
        test_number = yes_bound + max(1, (no_bound - yes_bound)//2)
        return search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop)

def best_topgraph_z_ini_mapping(L1, C, G, anchor, stop):
    ''' Test a more efficient way for scanning the gates in C than 'y' by using search_bipartite_top_z'''      
    g = graph_of_circuit(C)
    test = is_embeddable(g, G, anchor, stop)
    if test[0]:
        #print('The graph of the circuit is embeddable in G')
        return g, test[1]
    #print('The graph of the circuit is very likely NOT embeddable in G')    
    L_temp = [] #the index list of first cnot for each edge 
    for edge in g.edges():
        p, q = edge[0], edge[1]
        s = min([k for k in L1 if set(C[k]) == {p,q}]) 
        L_temp.append(s)

    L_temp.sort()
    yes_bound = 0
    no_bound = len(L_temp)
    test_number = no_bound//2
    exact_bound = search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop)
            
    g = nx.Graph()
    # add the first edge into g
    for s in L_temp[:exact_bound+1]:
        g.add_edge(C[s][0], C[s][1])
    test = is_embeddable(g, G, anchor, stop)
    if not test[0]: raise Exception('Check why the subgraph is not embeddable!')
    return g, test[1]

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# the topgraph initial mapping
def topgraph(C, G, anchor, stop):
    ''' Return the topgraph initial mapping

    Args:
        C (list): the input circuit
        G (graph): the architecture graph
    Returns:
        tau (list): the topgraph initial mapping
    '''    
    L = list(range(len(C))) # o in {o,x,y,z}
    tau_dict = best_topgraph_z_ini_mapping(L, C, G, anchor, stop)[1]
    
    return tau_dict

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

# Weighted SUBGRAPH
def best_wtg_o_ini_mapping(C, G, anchor, stop): #'o' for original
    ''' Return a graph g which is isomorphic to a subgraph of G
            while maximizing the number of CNOTs in C that correspond to edges in g
        Method: sort the edges according to their weights (the number of CNOTs in C corresponding to each edge);
                construct a graph by starting with the edge with the largest weight; then consider the edge with the second large weight, ...
                if in any step the graph is not isomorphic to a subgraph of G, skip this edge and consider the next till all edges are considered.
    
    Args:
        C (list): the input circuit
        G (graph): the architecture graph
        
    Returns:
        g (graph)
        map (dict)
    '''    
    g_of_c = graph_of_circuit(C)
    test = is_embeddable(g_of_c, G, anchor, stop)
    if test[0]:
        #print('The graph of the circuit is embeddable in G')
        return g_of_c, test[1]
    
    edge_wgt_list = list([C.count([e[0],e[1]]) + C.count([e[1],e[0]]), e] for e in g_of_c.edges())
    edge_wgt_list.sort(key=lambda t: t[0], reverse=True) # q[0] weight, q[1] edge
    
    '''Sort the edges reversely according to their weights''' 
    EdgeList = list(item[1] for item in edge_wgt_list)    
    #edge_num = len(EdgeList)
    
    '''We search backward, remove the first edge that makes g not embeddable, 
            and continue till all edges are evaluated in sequence. '''
            
    #Hard_Edge_index = 0 # the index of the first hard edge
    g = nx.Graph()
    result = dict()
    # add the first edge into g
    edge = EdgeList[0]
    g.add_edge(edge[0], edge[1])
    
    # rp = 0
    # for rp in range(edge_num):
    #     # h is the index of the last edge that can be added into g
    #     h = Hard_Edge_index
    #     if h == edge_num - 1: 
    #         return g, result
        
    #     EdgeList_temp = EdgeList[h+1:edge_num]
    #     for edge in EdgeList_temp:           
    #         g.add_edge(edge[0], edge[1])           
    #         i = EdgeList.index(edge)            
    #         # find the largest i such that the first i-1 edges are embeddable
    #         test = is_embeddable(g, G, anchor, stop)
    #         if not test[0]:
    #             Hard_Edge_index = i
    #             g.remove_edge(edge[0], edge[1])
    #             break
    #         result = test[1]
    #         if i == edge_num- 1:                
    #             return g, result
    # return g, result
    

    #EdgeList_temp = EdgeList[:]
    for edge in EdgeList:           
        g.add_edge(edge[0], edge[1])           
        test = is_embeddable(g, G, anchor, stop)
        if not test[0]:
            g.remove_edge(edge[0], edge[1])
            if nx.degree(g, edge[0]) == 0: g.remove_node(edge[0])
            if nx.degree(g, edge[1]) == 0: g.remove_node(edge[1])
        else:
            result = test[1]
    return g, result

# the weighted subgraph initial mapping
def wgtgraph(C, G, anchor, stop):
    ''' Return the weighted subgraph initial mapping

    Args:
        C (list): the input circuit
        G (graph): the architecture graph
    Returns:
        tau (list): the weighted subgraph initial mapping
    '''
    result = best_wtg_o_ini_mapping(C, G, anchor, stop)[1] # o in {o,x}
    return result 