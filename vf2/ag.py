#the architecture graphs of quantum devices

import networkx as nx

# IBM Q Tokyo (Q20) 
def q20():
    g = nx.Graph()
    g.add_nodes_from(list(range(20)))
    for i in range(0,4):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)

    for i in range(5,9):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(10,14):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(15,19):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(0,15):
        g.add_edge(i,i+5)
        g.add_edge(i+5,i)

    for i in [1,3,5,7,11,13]:
        g.add_edge(i,i+6)
        g.add_edge(i+6,i)

    for i in [2,4,6,8,12,14]:
        g.add_edge(i,i+4)
        g.add_edge(i+4,i)
    return g

#Grid AGs 
def qgrid(m,n):
    g = nx.Graph()
    g.add_nodes_from(list(range(0,m*n-1)))
    for i in range(0,m):
        for j in range(0,n):
            if i < m-1: g.add_edge(j*m+i, j*m+i+1)
            if j < n-1: g.add_edge(j*m+i, (j+1)*m+i)
    return g

# IBM Q Rochester (53 qubit)
def rochester():
    g = nx.Graph()
    g.add_nodes_from(list(range(0,53)))
    I_1= list(range(4)) + list(range(7,15)) +\
        list(range(19,27)) + list(range(30,38)) +\
            list(range(42,50))
            
    for i in I_1:
        g.add_edge(i,i+1)
    E = [(0,5),(5,9),(4,6),(6,13),(7,16),(16,19),\
         (11,17),(17,23),(15,18),(18,27),(21,28),(28,32),\
             (25,29),(29,36),(30,39),(39,42),(34,40),(40,46),\
                 (38,41),(41,50),(44,51),(48,52)]
    g.add_edges_from(E)
    return g

#Google Sycamore (53 qubits. Note this is relabelled so that the node set is continuous!)
def sycamore():
    g = nx.Graph()
    g.add_nodes_from(list(range(0,54))) 
    I = list(range(6,12))+list(range(18,24))+list(range(30,36))+\
        list(range(42,48))
    for i in I:
        for j in g.nodes():
            if j in I: continue
            if i-j in [5,6] or j-i in [6,7]:
                g.add_edge(i,j)
    g.remove_node(3)
    sorted(g)
    if 3 in g.nodes(): raise Exception('node error')
    mapping = dict()
    for n in g.nodes():
        if n < 3:
            mapping[n] = n
        else:
            mapping[n] = n - 1
            
    h = nx.relabel_nodes(g, mapping)
    sorted(h)
    return h

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#_________________TEST_____________________________________________________________#
if __name__=='__main__':
    # AG = sycamore()
    AG = qgrid(2,3)
    print(AG.edges())
    # print(nx.adjacency_matrix(AG,weight=None))
