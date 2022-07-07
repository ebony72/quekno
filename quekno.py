'''AG stands for the architecture graph of a quantum device '''
'''node_set of the AG should be continuous and start from 0, i.e., node_set = [0,k]'''
'''Permutations (except the 1st one) are induced by ≤2 swaps or by a set of parallel swaps'''

import networkx as nx
import random
import math
import copy
import csv
from qiskit import QuantumCircuit
from qctCircuit import qctCircuit
from utils_now import get_nodelist, is_parallel, permuted_sublist, is_embeddable,\
    action_to_perm, swap_circ, inverse_permutation, permuted_circuit,\
        permutation_affected_gates, map2tau

metapath = 'metaX/'
path = 'benchmarkX/'
    
class glinkChain:
    '''(1) generate a list of k sublists, 
       (2) generate a link (permutation) between any two consecutive sublists 
            such that the permutation is an embedding-destroyer,
       (3) generate an AG-circuit from each sublist and join them using the permutations 
    '''
    '''#_1qbg+2*#_2qubg = (qbg_ratio+2)*#_2qubg <= #device_node * 50 '''
    def __init__(self,
                 ag: nx.Graph,
                 subgraph_size: str, #'small', 'large'
                 cost_type: str, #'gate_opt1', 'gate_opt2' (2 swaps in each action) or 'depth_opt'
                 cost: int, #0, 5, 10, 15, 20, 25, 30, 35, 40 SWAPs
                 qbg_ratio, #QSE (0.51,0.4)=>2.55 TFL (0.27,0.36)=>1.5
                 circno: int):
        # internalizing initial parameters
        self.device_layout = ag
        self.connection_list = nx.edges(ag) 
        self.node_set = nx.nodes(self.device_layout) 
        self.subgraph_size = subgraph_size # 'large' or 'small'
        self.cost_type = cost_type # 'depth_opt', 'gate_opt1', 'gate_opt2'
        self.cost = cost #an upper bound of SWAP gates/depth required for QCT  
        self.qbg_ratio = qbg_ratio #QSE (0.51,0.4)=>2.55 TFL (0.27,0.36)=>1.5
        self.circno = circno # circuit_number
        self.exact = True
        self.name = "{}QBT_{}_{}_{}_{}_no.{}".format(len(self.node_set),subgraph_size,cost_type, cost, qbg_ratio, circno)
        self.num_qubit = None
               
        self.chain = {} # a dict of sections, each containing a set of AG-edges
        self.perm = {} # a dict of permutations #! Note these are not multiplied!
        self.action = {} # inserted swaps 
        self.destroy_2qbg = {} # 2qbg that are modified by the action
        self.permcirc = {} # a dict of permuted circuits 
        self.agcirc = {} # a dict of ag circuits
        
        self.depth_before, self.depth_after = 0, 0
        self.quekno_circ = [] #the generated quekno circuit
        self.transformed_circ = [] #the transformed circuit, i.e., an AG-circuit
        self.opt_transformed_circ = [] # a better transformed circuit
        self.quekno_sol = dict() #for each chain key, associate a pair (l2p_mapping,action)
        
        #procedure
        self.generate_glink_chain() 
        
        self.n_sec = len(self.chain) #number of sections        
        if self.chain.keys() != self.perm.keys():
            raise ValueError(self.chain, self.perm)

        #procedure        
        self.generate_quekno_circuit()

        self.get_quekno_sol()  #get the quekno solution of the quekno circuit  

        '''Determine the depths of the quekno circuit and its ideal transformation'''
        self.depth_before, self.depth_after = self.get_depth()
        
        '''verify if the quekno circuit has a solution with the specified cost'''
        self.verify() 
        
        # self.glc_print()

    
    def generate_random_sublist(self, k):
        '''generate a sublist with k differnt edges (from an AG)'''
        
        if k <= 0 and type(k) != int: 
            raise ValueError(k,type(k))
        num_added = 0
        sublist = []
        source_list = list(self.connection_list)
        while num_added < k:
            if not source_list:
                return sublist
            newitem = random.choice(source_list)
            source_list.remove(newitem)
            sublist.append(newitem)
            num_added += 1
        if len(source_list) >= len(self.connection_list): 
            raise ValueError(len(source_list), num_added, k)
        return sublist
    
    def single_edge_list(self, node_sublist):
        '''All single SWAPs relevant to the current sublist, here a SWAP is represented by an edge with form (p,q)'''

        temp_list = list(edge for edge in self.connection_list if edge[0] in node_sublist or edge[1] in node_sublist)
        return temp_list        
        
    def double_edge_list(self, node_sublist):
        '''All double SWAPs with form (x,y), where x,y are edges at least one of which are relevant to the current sublist'''
    
        sel = self.single_edge_list(node_sublist)[:]
        temp_list1 = list((x,y) for x in sel for y in self.connection_list if set(x) != set(y))
        temp_list2 = list((y,x) for x in sel for y in self.connection_list if set(x) != set(y))
        return temp_list1+temp_list2
    
    def generate_random_parallel_edge_list(self): #depth-opt
        '''generate a random parallel edge list'''

        EDGE  = list(self.connection_list)
        edge = random.choice(EDGE)
        action = [edge]
        EDGE.append('nonmax')
        hist = []
        while EDGE:
            node_action = set(get_nodelist(action))
            EDGE = list(edge for edge in EDGE if not set(edge)&node_action)
            edge = random.choice(EDGE)
            if edge in action or (edge[1],edge[0]) in action: continue
            if edge == 'nonmax': break 
            action.append(edge)
    
        #assertion 
        if not is_parallel(action): 
            print('**',EDGE)
            print('**',hist)
            raise ValueError('action is not parallel', action)
        return action
    
    #the depth of each subcircuit needs to be aligned 
    def circuit_align(self, sublist, C):  
        QC = qctCircuit(self.device_layout, C)
        d = QC.depth
        QDepth = QC.QubitDepth #dict
        NonMax = [q for q in QC.qubits if QC.QubitDepth[q] < d-1]
            
        while NonMax: 
            sublist_now = [edge for edge in sublist if set(edge) <= set(NonMax)]
    
            if not sublist_now:
                while NonMax:
                    p = random.choice(NonMax)
                    C.append((p,))
                    QDepth[p] += 1
                    if QDepth[p] == d-1: 
                        NonMax.remove(p)
                return C
                        
            gate_type = random.choice([1,2])
            if gate_type == 2: 
                gate = random.choice(sublist_now)
                C.append(gate)
                p,q = gate
                d_now = max(QDepth[p],QDepth[q])
                QDepth[p] = d_now+1
                QDepth[q] = d_now+1
                if QDepth[p] == d-1: 
                    NonMax.remove(p)
                    NonMax.remove(q)
            else:
                p = random.choice(NonMax)
                C.append((p,))
                QDepth[p] += 1
                if QDepth[p] == d-1: 
                    NonMax.remove(p)                
        return C 
    
    def generate_random_circuit(self, sublist, destroy_2qbg):
        '''generate a random circuit s.t. all cnots are from sublist and 1-qubit gates obey qbg_ratio'''
        
        node_sublist = get_nodelist(sublist)
        
        '''Put gates in destroy_2qbg front of the output circuit'''
        C = destroy_2qbg[:]
        count_1qbg_a = math.ceil(len(C) * self.qbg_ratio)
        count_1qbg_a_now = 0
        while count_1qbg_a_now < count_1qbg_a:
            gate = (random.choice(node_sublist),) #1qbg has form (p,)
            C.append(gate)
            count_1qbg_a_now += 1
        random.shuffle(C)
    
        '''Ensure every edge in sublist appears in the output circuit'''
        C_add = [gate for gate in sublist if gate not in destroy_2qbg]
        
        x = random.choice([0.05*i for i in range(1,5)]) 
        count_2qbg = math.ceil(len(sublist)*(1+x))
        count_2qbg_now = len(sublist)    
        while count_2qbg_now < count_2qbg:
            gate = random.choice(sublist)
            C_add.append(gate)
            count_2qbg_now += 1
            
        count_1qbg = math.ceil(count_2qbg * self.qbg_ratio)
        count_1qbg_now = count_1qbg_a
        while count_1qbg_now < count_1qbg:
            gate = (random.choice(node_sublist),)
            C_add.append(gate)
            count_1qbg_now += 1
            
        random.shuffle(C_add)
        C = C+C_add
        
        '''The last layer of C may involve only a small subset of the qubits. 
            Thus a fill-up process is required when the aim is to optimise depth''' 
        if self.cost_type == 'depth_opt':
            C = self.circuit_align(sublist,C)
        return C 
        
    #the first permutation acts as the inverse of the initial mapping and applies on the whole circuit        
    def generate_first_permutation(self): #bad sycamore node may be used if node_set is not continuous
        X = list(self.node_set)
        random.shuffle(X)
        pi = {i:X[i] for i in self.node_set}
        self.perm[0] = pi
        
    def destroy_embedding(self, perm, sublist1, sublist2):
        '''Check if perm is an embedding destroyer'''
        sublist = sublist1 + permuted_sublist(sublist2,perm)
        g = nx.Graph()
        g.add_edges_from(sublist)
        embeddable, _ = is_embeddable(g,self.device_layout,False,10)
        if embeddable :
            '''perm is NOT an embedding destroyer for the two sublists'''
            return False
        return True        
    
    def generate_link(self, sublist1, sublist2):
        '''generate a random permutation that destroys the embedding'''
        '''If fail then consider a different sublist2'''
        node_sublist1 = get_nodelist(sublist1)
        
        if self.cost_type == 'depth_opt':
            for  repeat in range(10): #repeat 10 times 
                action = self.generate_random_parallel_edge_list()
                perm = action_to_perm(self.node_set,action)
                if self.destroy_embedding(perm, sublist1, sublist2):
                    return 'success', action, perm
            return 'fail', action, perm 
        
        SWAP = self.single_edge_list(node_sublist1)
        ACT_PERM = []
        for edge in SWAP: #edge has form (p,q)
            perm = action_to_perm(self.node_set,(edge,))
            ACT_PERM.append([(edge,),perm])        

        if self.cost_type == 'gate_opt2': 
            SWAP2 = self.double_edge_list(node_sublist1)
            ACT_PERM2 = list([action, action_to_perm(self.node_set,action)] for action in SWAP2)
            ACT_PERM += ACT_PERM2
        
        random.shuffle(ACT_PERM)       

        for item in ACT_PERM:
            action, perm = item
            if self.destroy_embedding(perm, sublist1, sublist2):
                return 'success', action, perm            
        action, perm = random.choice(ACT_PERM)
        return 'fail', action, perm
    
    def random_generate_graph_edge_number(self):
        graph_size = math.ceil(len(self.node_set)/4) 
        if len(self.node_set) == 20:
            graph_size = 5 # tokyo has many connections
        if len(self.node_set) == 53:
            graph_size = 8
        if self.subgraph_size == 'small':
            return graph_size
        return 2*graph_size
                
    def generate_glink_chain(self):
       
        cost_now = 0
        step = 0
        subgraph_size = self.random_generate_graph_edge_number() # average size of the subgraphs
        k = math.ceil(random.gauss(subgraph_size,1))
        if k <= 0: k = 1
        
        self.chain[0] = self.generate_random_sublist(k)
        self.generate_first_permutation()
        self.action[0] = []
        self.destroy_2qbg[0] = []
        
        while cost_now < self.cost:
            step += 1
            '''Ensure the added link destroys the embedding'''
            Success = False
            repeat = 0 #try ten times if not succeed
            while not Success:  #Repeat until success
                repeat += 1
                k = math.ceil(random.gauss(subgraph_size,1))
                if k <= 0: k = 1
                self.chain[step] = self.generate_random_sublist(k)
                Suc, action, newperm = self.generate_link(self.chain[step-1],self.chain[step])
                
                #assertion
                if self.cost_type=='depth_opt' and not is_parallel(action):
                    print('The action is not parallel', step, Suc, action, repeat)
                    continue
                
                if Suc == 'success': break
                if repeat == 10: 
                    '''//todo: introduce customised techniques, 
                        such as generate a triangle or a cycle 5 to break the isomorphism!'''
                    break
            
            if Suc == 'fail': 
                self.exact = False
                print('Failed, and the cost is not exact',self.name,repeat,action)
            
            if self.cost_type == 'depth_opt':
                cost_now += 1
            else: #'gate_opt'
                cost_now += len(action)
            self.perm[step] = copy.copy(newperm)
            self.action[step] = copy.copy(action)
            self.destroy_2qbg[step] = permutation_affected_gates(self.perm[step], self.chain[step])

    def generate_quekno_circuit(self):
        Circ = dict()  #AG-circuit for each section      
        for i in self.chain.keys(): #each Circ[i] is an AG-circuit
            Circ[i] = self.generate_random_circuit(self.chain[i], self.destroy_2qbg[i])
            self.agcirc[i] = (Circ[i])[:] #physical circuit
            self.permcirc[i] = (Circ[i])[:] #logical circuit

        C = []
        for i in self.chain.keys():
            idx = self.n_sec -i -1 #permute backward
            C = permuted_circuit(Circ[idx]+C, self.perm[idx])
            
            # permute the sections 
            for j in self.chain.keys():
                if j < idx: continue
                self.permcirc[j] = permuted_circuit(self.permcirc[j], self.perm[idx])
           
        self.quekno_circ = C[:]

        
    def get_quekno_sol(self):
        #a qct_sol is a dict {0:l2pmapping,1:action1, 2:action2, ...}
        quekno_sol = dict()
        quekno_sol[0] = inverse_permutation(self.perm[0]) #a l2p_mapping from LC to PC

        for key in self.perm:
            if key == 0: continue
            # '''Achtung! The action needs not reverse!''' 
            quekno_sol[key] = self.action[key][:]
            
        self.quekno_sol = quekno_sol   

    def get_depth(self):
        '''get the depth of the transformed circuit'''
        '''//todo: Note this is sometimes not ideal: we may adjust the position before and after the inserted swaps!'''
        TC = copy.copy(self.agcirc[0]) #transformed circuit obtained by inserting the inverse actions
        for i in self.chain.keys():
            if i == 0: continue
            '''action was used to permute AG-subcircuit, thus reverse is nece. in order to get the transformed circ'''
            TC += swap_circ(self.action[i]) 
            TC += self.agcirc[i]
        self.transformed_circ = TC[:]

        QC = qctCircuit(self.device_layout,self.quekno_circ)
        self.num_qubit = QC.num_qubit
        
        opt_transformed_circ1 = QC.qct_output(self.quekno_sol, swapop=False)        
        self.opt_transformed_circ = QC.qct_output(self.quekno_sol, swapop=True)
        
        TQC0 = qctCircuit(self.device_layout,self.transformed_circ)
        TQC1 = qctCircuit(self.device_layout,opt_transformed_circ1)
        TQC2 = qctCircuit(self.device_layout,self.opt_transformed_circ)
        
        # print(TQC0.depth, TQC1.depth, TQC2.depth)
        
        if TQC0.depth < TQC2.depth: 
            print('The optimised output circuit has larger depth!', QC.depth, TQC0.depth, TQC1.depth, TQC2.depth, self.cost)
        return QC.depth, min(TQC0.depth, TQC1.depth, TQC2.depth)    
 
        
    def verify(self):
        
        '''verify if the quekno circuit is the collection of the permuted circuits''' 
        test_C = []
        for i in range(self.n_sec):
            test_C += self.permcirc[i]
        if self.quekno_circ !=  test_C:
            raise ValueError('They are not equal', len(test_C),len(self.quekno_circ))
                
        '''Verify if the output circuit has a solution with cost close to self.cost'''
        if self.cost_type == 'depth_opt':
            print(self.depth_after, self.depth_before, self.cost, self.depth_after-self.depth_before-self.cost*3)
            if self.depth_after > self.depth_before+self.cost*3:
                print('The quekno circuit has real depth cost larger than', self.cost)
                self.exact = False
        
        '''Verify if the quekno_sol is a qct solution of the quekno circuit'''
        QC = qctCircuit(self.device_layout,self.quekno_circ)        
        QC.qct_verify(self.quekno_sol)
        
        '''Check if C is an AG-circuit'''
        C = copy.copy(self.opt_transformed_circ)
        for gate in C:
            if len(gate)==1: continue
            p,q = gate
            if (p,q) not in self.connection_list:
                print(p,q, (p,q) in nx.edges(self.device_layout), (p,q) in self.connection_list)
                raise ValueError('The transformed circuit is not executable!', gate)
                

    def glc_print(self):
        print('Quekao circuit :: \n', self.quekno_circ)
        print('Transformed circuit :: \n', self.transformed_circ)
        print('Opt-Transed circuit :: \n', self.opt_transformed_circ)
        print('Quekao chain :: \n', self.chain)
        print('agcirc :: \n', self.agcirc)
        print('permcirc :: \n', self.permcirc)
        print('Quekao permutation :: \n', self.perm)        
        print('Quekao action :: \n' , self.action)
        print('QCT solution :: \n', self.quekno_sol)
        
        for key in self.chain:
            print('QCT solution ::', key, self.quekno_sol[key])
        for key in self.chain:
            print('Quekao action ::',key, self.action[key])
            print('Quekao permutation ::', key, self.perm[key])        
            print('Quekao chain ::',  key, self.chain[key])
            print('permcirc ::', key, self.permcirc[key])        
            
    def save_quekno_sol(self): 
        '''The writer transcripts {3: 0, 1: 1, 5: 2, 4: 3, 0: 4, 2: 5} to [3, 1, 5, 4, 0, 2]'''
        '''i.e., the l2p_mapping is inversed and should change back from the csvFile''' 
        with open(metapath+self.name + "_solution.csv", 'w') as csvFile:
            writer = csv.writer(csvFile)
            for i in range(self.n_sec):
                if i==0: 
                    '''inimap is a complete l2p_mapping, tau is a list'''
                    tau = map2tau(self.quekno_sol[0],self.node_set())
                    writer.writerow(tau)
                else:
                    writer.writerow(self.quekno_sol[i])
        csvFile.close()
        
    def save_circ_as_qasm(self,circ,type):        
        qc = QuantumCircuit(len(self.node_set)) #qiskit quantum circuit
        for gate in circ:
            if len(gate) == 1:
                p = gate[0]
                qc.h(p)
            else:
                p,q = gate
                qc.cx(p,q)
        if type==0:  #quekno circ
            name = self.name
        elif type == 1:
            name = '[transformed]_' + self.name 
        else: #type == 2
            name = '[opt_transed]_' + self.name
        qc.qasm(formatted=False, filename = path + name + '.qasm')
                    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if __name__=='__main__':
    import vf2.ag# architecture graph
    '''Achtung! The sycamore in ag is relabelled and different from that used in QUEKO!'''
    from qctCircuit import save_record
    # AG = vf2.ag.qgrid(2,3)
    # AG = vf2.ag.q20()
    # AG = vf2.ag.sycamore()
    AG = vf2.ag.rochester()
    EG = list(AG.edges())
    # print(AG.edges())
    
    def costset(cost_type):
        if cost_type == 'depth_opt':
            cost_set = {1,2,3,4,5,10}
        else:
            cost_set = {0,1,2,3,4,5,10,15,20,25} 
        return cost_set
    def parameter_combination():
        '''As tokyo is highly connected, we consider only large subgraphs'''
        '''For gate optimality, only cx gates are concerned, so qbg_ratio is irrelevant'''
        Para_Com = []
        for subgraph_size in {'small','large'}: 
            if len(nx.nodes(AG)) <= 30 and subgraph_size == 'small': continue
            for cost_type in {'gate_opt1', 'gate_opt2', 'depth_opt'}:
                cost_set = costset(cost_type)
                for cost in cost_set:
                    for qbg_ratio in {1.5, 2.55}:
                        if cost_type != 'depth_opt' and qbg_ratio == 2.55: continue
                        Para_Com.append([subgraph_size,cost_type,cost,qbg_ratio])        
        return Para_Com

    '''For ag=tokyo and gate_opt, we have 1*2*10*1*10=200 benchmark circuits'''
    '''For ag=tokyo and depth_opt, we have 1*1*6*2*10=120 benchmark circuits'''
    '''For ag=sycamore/rochester and gate_opt, we have 2*2*10*1*10=400 benchmark circuits'''
    '''For ag=sycamore/rochester and depth_opt, we have 2*1*6*2*10=240 benchmark circuits'''
    
    record = []
    record.append(('name','cost_type','cost','num_qbt', 'exact', 'gate_in', 'gate_out',\
                   'cx_in','cx_out','cx_ratio','depth_in','depth_out','depth_ratio'))
    
    Para_Com = parameter_combination()
    for combintuple in Para_Com:
        subgraph_size, cost_type, cost, qbg_ratio = combintuple
        
        '''Testing! generate two benchmarks, one with small subgraph, the other large'''
        if cost_type != 'gate_opt2' or cost != 10 or qbg_ratio != 1.5: continue
        for circno in range(0,1): 
        # for circno in range(0,10): 
            gc = glinkChain(AG,subgraph_size, cost_type, cost, qbg_ratio, circno)
            
            # gc.glc_print()
            '''If you need to output the benchmarks as qasm files ...'''
            gc.save_quekno_sol()
            gc.save_circ_as_qasm(gc.quekno_circ,0)        
            # gc.save_circ_as_qasm(gc.transformed_circ,1)        
            # gc.save_circ_as_qasm(gc.opt_transformed_circ,2)    

            cxcirc = [gate for gate in gc.quekno_circ if len(gate)==2]
            out_cxcirc = [gate for gate in gc.transformed_circ if len(gate)==2]
            record.append((gc.name,gc.cost_type,gc.cost,gc.num_qubit, gc.exact,\
                           len(gc.quekno_circ), len(gc.transformed_circ),\
                               len(cxcirc), len(out_cxcirc), round(len(out_cxcirc)/len(cxcirc),4),
                                   gc.depth_before, gc.depth_after, round(gc.depth_after/gc.depth_before,4)))

    name =  '53QBT_'+cost_type+'_Rochester_'+ 'quekno_benchmark_info'
    save_record(metapath, name, record)

