#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  7 07:27:59 2022

@author: sanjiangli
"""

'''output the basic info of qeukao benchmarks// this is not required as we have generated 
    benchmark info in quekno.py. It can serve as a double checking purpose or when you want to 
    extract the qct solution. 
'''


import os, csv, re
import vf2.ag
from utils_now import CreateCircuitFromQASM, ReducedCircuit, save_record
from qctCircuit import qctCircuit


'''A quekno solution is a dictionary where the value for key=0 is the initial mapping;
    the other values are actions [edge1, edge2, ...] ''' 
def get_qct_sol(path, filename): #from ZXZ
    '''read the quekno soluation from the csv file'''
    '''Note in the csv file, a l2p_mapping like {3: 0, 1: 1, 5: 2, 4: 3, 0: 4, 2: 5}\
        was transcripted to [3, 1, 5, 4, 0, 2] and should be reversed'''
    
    quekno_sol = dict()
    with open(path+filename, newline='') as csvfile:
        queknoreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in queknoreader:        
            idx = queknoreader.line_num - 1            
            if idx == 0:
                row = row[0]
                row += ','
                datas = re.findall("([0-9]+),", row)
                row_convert = [int(data) for data in datas]
                quekno_sol[0] = dict()
                for i in row_convert:
                    quekno_sol[0][row_convert[i]] = i
                # print(idx,row_convert)
            else:
                row_convert_str = ''
                for r in row: row_convert_str += r
                datas = re.findall("\(([0-9]+,[0-9]+)\)", row_convert_str)
                datas1 = re.findall("\(([0-9]+),", row_convert_str)
                datas2 = re.findall(",([0-9]+)\)", row_convert_str)
                row_convert = []
                for data1, data2 in zip(datas1, datas2):
                    row_convert.append((int(data1), int(data2)))
                quekno_sol[idx] = row_convert
                # print(idx,row_convert)
    return quekno_sol

        
'''The mappings in the solutions are p2l mappings, so need to transfer to l2p mappings'''


# name = 'quekno_info_' + 'test_' 
# AG = vf2.ag.qgrid(2,3)
path = "./benchmarkX/"
metapath = './metaX/'

# name = 'quekno_info_depth_' + 'tokyo_test_' 
# AG = vf2.ag.q20()

AG = vf2.ag.rochester()
name = 'quekno_info_' + 'rochester_test_' 


# AG = vf2.ag.sycamore()
# name = 'quekno_info_gate_' + 'sycamore_test_' 


node_num = len(AG.nodes())
label = [x for x in range(node_num)]
        
record = []
record.append(('name','num_qbt','num_gate','num_cx','num_+cx','num_+real_cx','num_depth','num_+depth','num_+real_depth'))

sum_depth_in, sum_depth_out = 0, 0
for filename in os.listdir(path):
    if not filename.endswith('.qasm'): continue
    print(filename)
    
    '''Extract the circuit from qasm files'''
    circ_qasm = CreateCircuitFromQASM(filename, path)
    circ = ReducedCircuit(circ_qasm) #1&2qbg list
    

    QC = qctCircuit(AG, circ)
    qubit_num_circ = QC.num_qubit 
    cx_gates = QC.cxcirc
    
    sol_filename = filename[0:-5]+'_solution.csv' 
    quekno_sol = get_qct_sol(metapath, sol_filename)

    QC.qct_verify(quekno_sol)   
        
    '''output circuit'''
    output_circ = QC.qct_output(quekno_sol,swapop=False)        
    OQC = qctCircuit(AG, output_circ)
    
    '''optimised output circuit, which may have shallow depth but not always so''' 
    output_circ_op = QC.qct_output(quekno_sol,swapop=True)        
    OQC_op = qctCircuit(AG, output_circ_op)
    
    num_cx = QC.num_cx
    print(filename, num_cx, QC.depth,OQC.depth,OQC_op.depth)
    num_cxplus = OQC.num_cx - num_cx 
      
    output_depth = min(OQC.depth,OQC_op.depth)
    content = (filename, qubit_num_circ, QC.num_gate, num_cx, num_cxplus,\
                round(OQC.num_cx/num_cx,4), QC.depth, output_depth, round(output_depth/QC.depth,4))
    print(content)
    record.append(content)
    # save_record(name, record)
    
    sum_depth_in += QC.depth
    sum_depth_out += output_depth
print(sum_depth_in, sum_depth_out, round(sum_depth_out/sum_depth_in, 4))
        
    