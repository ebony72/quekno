#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 09:53:08 2022

@author: sanjiangli
"""
from vf2.vfs import Vf
import networkx as nx
import csv
import copy
from qiskit import QuantumCircuit 

def save_result(path, name, content):#//todo
    name = str(name)
    content = str(content)
    file = open(path + '_' + name, mode = 'a')
    file.write(content)
    file.write('\n')
    file.close() 

def save_record(path, name, record):
    with open(path + '_' + name + ".csv", 'w') as csvFile:
        writer = csv.writer(csvFile)
        for item in record:
            writer.writerow(item)
    csvFile.close()
    
spl = nx.shortest_path_length
def SPL(g):
    if not nx.is_connected(g): raise Exception ('g is not connected!')
        # '''g is disconnected! We consider its largest connected component instead!'''
        # largest_cc = max(nx.connected_components(g), key=len)
        # g = g.subgraph(largest_cc)
    spl_dic = dict() 
    V = list(g.nodes())
    V.sort()
    for p in V:
        for q in V:
            if (q,p) in spl_dic:
                d = spl_dic[(q,p)]
            else:
                d = nx.shortest_path_length(g,p,q)
            spl_dic[(p,q)] = d
    return spl_dic

#assertion
def is_parallel(action): #action is a seq. of edges/swaps
    if len(action)<2: return True
    for edge in action:
        if action.count(edge) > 1: return False
    for edge1 in action:
        for edge2 in action:
            if edge1 == edge2: continue
            if set(edge1)&set(edge2): 
                return False
    return True
    
def is_embeddable(g, H, anchor, stop):
    '''//todo: could make this customised to H (the architecture graph)'''
    '''check if a small graph g is embeddable in a large H, anchor is bool
        g, H (Graph)
        anchor (bool): whether or not mapping anchor of g to that of H
        stop (float): time limit for vf2
    '''
    vf2 = Vf()
    result = {} 
    # if anchor: result[hub(g)] = hub(H)
    result = vf2.dfsMatch(g, H, result, stop)
    lng = len(nx.nodes(g))
    if len(result) == lng:
        return True, result   
    return False, result

def get_nodelist(sublist): # nodes in a list of cx gates
    node_list = list(x[0] for x in sublist) 
    node_list += list(x[1] for x in sublist if len(x)==2)
    node_list = list(set(node_list))
    return node_list

def reverse_action(action):
    reverse = action[:]
    if len(action) == 2:
        if len(action[0]) != 2: raise ValueError('!!', action)
        edge1, edge2 = action
        if set(edge1) & set(edge2):
            reverse = (edge2,edge1)
    return reverse    
    
def swap_circ(action):
    '''Produce the swap circuit of the reverse of action'''
    '''Here we only consider single_edge, double_edge and parallel_edge actions'''

    swap_circ = []
    for edge in action:
        p, q = edge
        swap_circ += ((p,q),(q,p),(p,q))
    return swap_circ

def find_reverse_mapping(mapping, qubit_num): #sl: revised@220503
	reverse_mapping = dict()
	for l_qubit in mapping:
		p_qubit = mapping[l_qubit]
		reverse_mapping[p_qubit] = l_qubit
	return reverse_mapping

#not used yet
def inverse_permutation(perm):
    newperm = copy.deepcopy(perm)
    if type(perm)==list:
        for i in range(len(perm)):
            newperm[i] = perm[i]
    if type(perm) == int:
        raise TypeError(perm)
    inv_perm = dict()
    for key in perm:
        val = perm[key]
        inv_perm[val] = key
    return inv_perm     

def action_to_perm(node_set,action): #depth-opt
    '''The permutation induced by a swap action (edge1,edge2,...)'''
    if not action: raise ValueError('action should not be empty!', action)
    perm = {i:i for i in node_set}
    for i in range(0,len(action)):
        u,v = action[i]
        perm[u],perm[v] = perm[v],perm[u]   #swap the value of u, v
    return perm

def permuted_sublist(sublist, perm): #0605 Here perm is a p2l_mapping
    '''Permutate a sublist with a permutation perm'''
    temp_list = []
    for edge in sublist:
        p,q = edge
        u,v = perm[p],perm[q]
        temp_list.append((u,v))
    return temp_list

def get_gate_idx_list_on_qubit(node_list,circ): #circ 1&2qbg list
    gates_on_qubit = dict()
    for p in node_list:
        gates_on_qubit[p] = []
    for idx in range(len(circ)):
        gate = circ[idx]
        if len(gate) == 2:
            p,q = gate
            gates_on_qubit[p].append(idx)
            gates_on_qubit[q].append(idx)
        else:
            p = gate[0]
            gates_on_qubit[p].append(idx)
    return gates_on_qubit
                    
def permuted_circuit(circ, perm): #0605 perm is a p2l_mapping
    '''Permutate a circuit with a permutation, where a gate in circ has form 
        either (p,) or (p,q)'''
    temp_circ = []
    for gate in circ:
        if len(gate) == 1:
            p = gate[0]
            u = perm[p]
            temp_circ.append((u,))
        else: #len(gate)==2
            p,q = gate
            u,v = perm[p],perm[q]
            temp_circ.append((u,v))
    return temp_circ

def permutation_affected_gates(perm, circ):
    '''collection of cx gates changed by perm'''
    C_aff = []
    for gate in circ:
        if len(gate) != 2: continue
        p,q = gate
        if {p,q} == {perm[p],perm[q]}: continue
        C_aff.append(gate)
    return C_aff


def CreateCircuitFromQASM(file, path):
    QASM_file = open(path + file, 'r')
    iter_f = iter(QASM_file)
    QASM = ''
    for line in iter_f: 
        QASM = QASM + line
    #print(QASM)
    cir = QuantumCircuit.from_qasm_str(QASM)
    QASM_file.close    
    return cir

def CreateQASMFromQuekaoCircuit(k, quekno_circ):
    qc = QuantumCircuit(k)
    for gate in quekno_circ:
        if len(gate) == 1:
            p = gate[0]
            qc.h(p)
        p,q = gate
        qc.cx(p,q)
    return qc

def ReducedCircuit(circ_qasm):
    '''Return Reduced Circuit containing only [name, [p,q]], e.g., ['cx', [0,2]] '''
    C = []
    for gate in circ_qasm:
        # if gate[0].name != 'cx': continue
        qubits = [q.index for q in gate[1]]
        C.append(qubits)
    # self.circuit = C[:]
    return C

def tau2map(tau):
    ''' Return the l2p_mapping: Q -> V of tau: V -> Q '''
    dic1 = dict()
    for q in tau:
        if q != -1: #i.e., q has been assigned the phy. qubit u = tau.index(q)
            dic1[q] = tau.index(q)
    return dic1

def map2tau(dic, V):
    '''Return the tau: V -> Q of a (l2p_mapping) dic: Q -> V '''
    tau = [-1]*len(V)
    for q in dic:
        tau[dic[q]] = q
    return tau

def swap(l2p_mapping, u, v, EG):

    if (u, v) not in EG: 
        raise ValueError('invalid edge',u,v)
    # print(l2p_mapping)
    p2l_mapping = inverse_permutation(l2p_mapping)

    if u in p2l_mapping and v in p2l_mapping:
        p2l_mapping[u], p2l_mapping[v] = p2l_mapping[v], p2l_mapping[u]
    elif u in p2l_mapping and v not in p2l_mapping:
        p2l_mapping[v] = p2l_mapping[u] 
        del p2l_mapping[u]
    elif u not in p2l_mapping and v in p2l_mapping:
        p2l_mapping[u] = p2l_mapping[v] 
        del p2l_mapping[v]
        
    l2p_mapping = inverse_permutation(p2l_mapping)
    return l2p_mapping     

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if __name__=='__main__':
    node_set = {0,1,2,3,4,5}
    action = ((0,2),(2,3))
    perm  = action_to_perm(node_set,action)
    print(perm)

