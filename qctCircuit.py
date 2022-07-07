#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 19:45:38 2022

@author: sanjiangli
"""
'''We have the following formats of circuits:
    - QASM
    - 1qbg&2qbg list: 1qbg with form (p,) and 2qbg with form (p,q)
        (for depth-opt qct, we need consider 1qbg and 2qbg together)
    - 2qbg list (p,q) or [p,q] (for gate-opt qct, we need only 2qbg list)
'''
import networkx as nx
import copy
from utils_now import *

class qctCircuit:
    '''ag (nx.Graph) : the architecture graph 
       circuit (list) : 1&2qb gate list'''
    def __init__(self, ag, circuit):
        self.device_layout = ag
        self.connection_list = nx.edges(ag)
        self.node_list = list(nx.nodes(ag))
        self.num_node = len(nx.nodes(ag))
        self.circuit = circuit #1qbg&2qbg form
        self.num_gate = len(self.circuit)
        self.IDX = [i for i in range(self.num_gate)] #index list
        self.cxcirc = self.ReducedCXCircuit()
        self.num_cx = len(self.cxcirc)
        self.qubits = self.qubit_in_circuit()
        self.num_qubit = len(self.qubits)
        self.output_circ = []
        
        self.Layer = dict() #partiton circ into layers
        self.QubitDepth = dict() #int
        self.depth = -1 #! if a qubit is not occupied, then its depth is -1
        self.cxgraph = nx.Graph()
        self.dpdgraph = nx.DiGraph() #dependency graph 
        self.Depth_Dic = dict()
        
        self.LayerPartition()      
        self.get_cxgraph() 
        self.get_dependency_graph(self.IDX, has_1qbt=True)
        
        for d in range(self.depth):
            for idx in self.Layer[d]:
                if self.Depth_Dic[idx] != d:
                    raise ValueError (d, idx, self.Layer[d], self.Depth_Dic[idx])
                    
    def ReducedCXCircuit(self):
        cxcirc = []
        for gate in self.circuit:
            if len(gate) != 2: continue
            cxcirc.append(gate)
        return cxcirc
    
    def qubit_in_circuit(self): 
        Qset = set()
        for gate in self.circuit:
            Qset.add(gate[0])
            if len(gate)==2:
                Qset.add(gate[1])
        return Qset

    def topgates(self, LD): 
        ''' Return the index list of top layer gates, LD is current gate idx list'''
        if not LD: return []
        LTG = []
        N = set() # the set of qubits of gates in LTG
        for i in LD:
            if N == self.qubits: return LTG                    
            gate = self.circuit[i]
            p = gate[0]
            if len(gate) == 1:
                if p in N: continue
                N.add(p)
                LTG.append(i)
            else:
                q = gate[1]
                if p in N and q in N: continue
                if p not in N and q not in N:
                    N.add(p)
                    N.add(q)
                    LTG.append(i)
                elif p not in N: #q in N
                    N.add(p)
                elif q not in N: #p in N
                    N.add(q)
            if N == self.qubits: return LTG        
        return LTG
        
    def LayerPartition(self):        
        Layer = dict()
        QubitDepth = dict() #the depth of qubit q
        for q in self.node_list:
            QubitDepth[q] = -1
            
        depth = 0
        LD = self.IDX[:]
        while LD:
            X = copy.deepcopy(self.topgates(LD))
            if not X: raise ValueError('topgates is empty', X)
            
            qubits_now = set()
            for i in X:
                gate = self.circuit[i]
                qubits_now = qubits_now.union(set(gate))

            for q in qubits_now:
                QubitDepth[q] = depth #update the depth of each qubit
                
            Layer[depth] = X
            for idx in X:
                LD.remove(idx)
                
            depth += 1
        self.Layer = Layer
        self.QubitDepth = QubitDepth
        self.depth = len(self.Layer)

    def get_cxgraph(self):   
        g = nx.Graph()
        g.add_nodes_from(self.qubits)
        for gate in self.circuit:
            if len(gate) != 2: continue
            g.add_edge(gate[0],gate[1])
        self.cxgraph = g          

    # compute the gate dependency graph of (a part of) the circuit
    def get_dependency_graph(self, LX, has_1qbt): 
        '''LX (list) is the index set of a subset of C, has_1qbt (bool)'''
        
        g = nx.DiGraph()
        g.add_nodes_from(LX)
        L1 = LX[:]
        Depth_dict = dict()
        depth = 0
        C = copy.deepcopy(self.cxcirc)
        if has_1qbt: C = copy.deepcopy(self.circuit)
        
        while L1:
            L0 = L1[:]
            TG = self.topgates(L0)
            for i in TG:
                L0.remove(i)
                Depth_dict[i] = depth                
            depth += 1
            if not L0: break                   
            for i in TG:
                p = C[i][0]
                for j in L0:
                    if p in C[j]:
                     g.add_edge(i,j) 
                     break
                if len(C[i]) == 2:
                    q = C[i][1]
                    for j in L0:
                        if q in C[j]:
                             g.add_edge(i,j) 
                             break
                         
            L1 = L0[:]
        self.dpdgraph = g 
        self.Depth_Dic = Depth_dict

    def entail(self, l2p_mapping, gate):
        ''' Check if an l2pmapping entails a gate
        Args:
            l2p_mapping (dict): a mapping from Q to V 
            gate (list): [p,q] represents a CNOT gate
        Returns:
            Boolean: True if l2p_mapping entails gate
        '''
        if len(gate) == 1: 
            p = gate[0]
            if p in l2p_mapping: return True
            return False
        p, q = gate[0], gate[1]
        if (p not in l2p_mapping) or (q not in l2p_mapping):
            return False
        u, v = l2p_mapping[p], l2p_mapping[q] 

        if (u,v) in self.connection_list: 
            return True    
        return False
        
    def greedy_solved_gates(self, l2p_mapping, LD):
        ''' Return the list of the indices of all gates solvable by the l2p_mapping
        Args:
            l2p_mapping (list): a mapping from Q to V, where V (Q) is the list of physical (logical) qubits
            LD (list): the index list of D, which represents the current logical circuit            
        Returns:
            LSTG (list): the list of indices of solvable gates in D
        '''    
        LDx = LD[:]
        LSTG = [] # the index set of all solved topgates 
        while True: # as long as there are gates that can be solved by this mapping
            if not LDx: return LSTG        
            LX = self.topgates(LDx) # topgates returns indices of topgates
            ldx = len(LDx)       
            # examine gates in LX one by one, check if they are entailed by l2p_mapping
            for i in LX:
                gate = self.circuit[i]
                if self.entail(l2p_mapping, gate):
                    LDx.remove(i)
                    LSTG.append(i) # C[i] is solved by l2p_mapping 
            if len(LDx) == ldx: # l2p_mapping does not solve any new topgate
                return LSTG

    def qct_verify(self, qct_sol):
        '''qct_sol (dict) '''
 
        #initial mapping
        l2p_mapping = copy.copy(qct_sol[0])  
        # print('initial mapping :: ', l2p_mapping)       

        LD = self.IDX[:]
        for key in qct_sol:
            if key > 0: 
                act = qct_sol[key]                        
                for edge in act:
                    # print('edge', edge)
                    u, v = edge
                    l2p_mapping = swap(l2p_mapping, u, v, self.connection_list)
                    # print('cur mapping :: ', l2p_mapping)       

            section = self.greedy_solved_gates(l2p_mapping, LD)            
            for idx in section:
                LD.remove(idx) 
                
        if LD: 
            print('Unsolved gates are in', [[idx,self.circuit[idx]] for idx in LD])
            raise ValueError('The solution is incorrect!')

    def qct_output(self, qct_sol, swapop=True):
        '''swapop (bool): ?optimisation is used to reduce depth of the ouput circuit''' 
        '''qct_sol (dict): 1st value = initial mapping and the rest actions'''
        
        temp_gate_idx_list = self.IDX[:]
        output_circ = []

        l2p_mapping = copy.copy(qct_sol[0]) #initial mapping       
        
        first_section = self.greedy_solved_gates(l2p_mapping, temp_gate_idx_list)
        for idx in first_section:
            temp_gate_idx_list.remove(idx)
        l_circ_sec_now = [self.circuit[idx] for idx in first_section]
        p_circ_sec_now = permuted_circuit(l_circ_sec_now, l2p_mapping)
        output_circ = p_circ_sec_now[:] #the first physical section
        
        if not swapop: #solve the gates greedily        
            for key in qct_sol:
                if key == 0: continue

                act = qct_sol[key]
                for edge in act:
                    # print('edge', edge)
                    u, v = edge
                    l2p_mapping = swap(l2p_mapping, u, v, self.connection_list)

                section = self.greedy_solved_gates(l2p_mapping, temp_gate_idx_list)
                for idx in section:
                    temp_gate_idx_list.remove(idx)

                action = copy.copy(qct_sol[key])
                '''The action = (edge1,edge2) needs not reverse as 
                    (i) in constructing quekao it is used to transform an AG-circuit 
                        first by swapping edge1, followed by swapping edge2
                    (ii) in verification it is used to transform a l2p_mapping 
                        first by swapping edge1, followed by swapping edge2.
                '''
                output_circ += swap_circ(action)
                
                l_circ_sec_now = [self.circuit[idx] for idx in section]    
                p_circ_sec_now = permuted_circuit(l_circ_sec_now, l2p_mapping)
                output_circ += p_circ_sec_now
            return output_circ
        
        '''The above construction is straightforward and could be depth-optimised.'''
        for key in qct_sol:
            if key == 0: continue
            '''move not aligned 1qb gates after the SWAP actions '''
            action = copy.copy(qct_sol[key])

            while len(action): 
                '''consider the first edge/swap in the action'''
                # print(action)
                v1,v2 = action[0]
                l2p_mapping = swap(l2p_mapping, v1, v2, self.connection_list)
                action = action[1:] #remove the swap we are examining 
                
                '''Consider a temporary qctCircuit'''     
                cur_output_circ = output_circ[:]
                tempOC = qctCircuit(self.device_layout,cur_output_circ)
                qubit_depth_dict = tempOC.QubitDepth
                
                '''Select not aligned 1qb gates'''
                gates_on_qubit = get_gate_idx_list_on_qubit(self.node_list, cur_output_circ) #dict
                '''Compare the progress at qubit v1 and v2'''
                diff  = qubit_depth_dict[v1] - qubit_depth_dict[v2] 
                
                if diff == 0: continue  
                if diff < 0:
                    u1,u2 = v1,v2
                else:
                    u1,u2 = v2,v1
                    
                Put_After_SWAP = []
                #push some gates on u2 forward (-->)
                gates_on_qubit[u2].reverse()
                examine_gate_list = gates_on_qubit[u2][0:abs(diff)]
                for idx in examine_gate_list:
                    gate = cur_output_circ[idx]
                    if u2 not in gate: 
                        raise ValueError('gates_on_qubit error', v2, gate,gates_on_qubit[v2])
                    if len(gate) == 2: #cannot push forward
                        break
                    output_circ.remove(gate)
                    Put_After_SWAP.append(idx)
                moved_gates = ((u1,),)*len(Put_After_SWAP) #we denote a 1qbg as (q,)
                
                output_circ += [(v1,v2),(v2,v1),(v1,v2)]
                output_circ += moved_gates
    
            cur_section = self.greedy_solved_gates(l2p_mapping, temp_gate_idx_list)
            # print(l2p_mapping)

            for idx in cur_section:
                # print('newly solved gate', idx, self.circuit[idx])
                temp_gate_idx_list.remove(idx)
                
            l_circ_sec_now = [self.circuit[idx] for idx in cur_section]
            p_circ_sec_now = permuted_circuit(l_circ_sec_now, l2p_mapping)
            output_circ += p_circ_sec_now
            
        if temp_gate_idx_list: 
            for key in qct_sol:
                print(key, qct_sol[key])
            raise ValueError('The solution is incorrect!')
        return output_circ
    
if __name__=='__main__':
    import ag # architecture graph
    '''Achtung! The sycamore in ag is relabelled and different from that used in QUEKO!'''
    AG = ag.qgrid(2,3)
    # AG = ag.q20()
    
    LC = [[1], [1, 0], [1], [2, 3], [2], [2], [3, 1], [3, 5]]
    qct_sol = {0: [5, 3, 0, 1, 4, 2], 1: ((2, 0),)}
    
    QC = qctCircuit(AG, LC)
    QC.qct_verify(qct_sol)
    output_circ = QC.qct_output(qct_sol,True)
    print(len(LC),len(output_circ))
    OQC = qctCircuit(AG, output_circ)
    print(QC.depth,OQC.depth)
    
    def save_quekno_circ_as_qasm(num_qubit, circ):        
        qc = QuantumCircuit(num_qubit) #qiskit quantum circuit
        for gate in circ:
            if len(gate) == 1:
                p = gate[0]
                qc.h(p)
            else:
                p,q = gate
                qc.cx(p,q)
        qc.qasm(formatted=False, filename = "benchmark/" + '_' +'test_output1.qasm')    

    save_quekno_circ_as_qasm(len(nx.nodes(AG)), output_circ)
    
    