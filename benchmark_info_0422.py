#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  7 07:27:59 2022

@author: sanjiangli
"""

'''output the basic info and qct-solutions of old qeukno benchmarks// only for 0422 version '''
'''We only have generated solutions for Tokyo and Rochester, Sycamore solutions were covered by those for Rochester.'''
''' The solution for `53QBT_gate_Rochester_large_opt1_10_1.5_no.7.qasm` was corrupted '''

import os, csv, re
import vf2.ag
from utils_now import CreateCircuitFromQASM, ReducedCircuit, save_record
from qctCircuit import qctCircuit


'''A quekao soluation is a dictionary where the value for key=0 is the initial mapping;
    the other values are actions [edge1, edge2, ...] ''' 
def get_qct_sol(path, filename): #from ZXZ
    '''read the quekao solution (0422 version) from the csv file'''
    '''The solution has name like `53QBT_depth_Rochester_large_opt_1_1.5_no.0_solution.csv` '''
    
    quekno_sol = dict()
    with open(path+filename, newline='') as csvfile:
        quekaoreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in quekaoreader:        
            idx = quekaoreader.line_num - 1            
            if idx == 0:
                row = row[0]
                row += ','
                datas = re.findall("([0-9]+),", row)
                row_convert = [int(data) for data in datas]
                quekno_sol[0] = dict()
                
                for i in row_convert:
                    quekno_sol[0][i] = row_convert[i] #for 0422 version
                # print(idx,row_convert)
                # print(quekno_sol[0])
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

AG = vf2.ag.rochester()
metapath = "./meta0422/rochester/"
name = 'quekao_info_58' + '_rochester_gate' 
path = "./benchmark0422/53Q_gate_Rochester/"
# name = 'quekao_info_58' + '_rochester_depth' 
# path = "./benchmark0422/53Q_depth_Rochester/"

# AG = ag.q20()
# metapath = "./meta0422/tokyo/"
# name = 'quekao_info_58' + '_tokyo_depth' 
# path = "./benchmark0422/20Q_depth_Tokyo/"
# name = 'quekao_info_58' + '_tokyo_gate' 
# path = "./benchmark0422/20Q_gate_Tokyo/"


node_num = len(AG.nodes())
label = [x for x in range(node_num)]
        
record = []
record.append(('name','num_qbt','num_gate','num_cx','num_+cx','num_+real_cx','num_depth','num_+depth','num_+real_depth'))

sum_depth_in, sum_depth_out = 0, 0
for filename in os.listdir(path):
    if not filename.endswith('.qasm'): continue
    if filename == '53QBT_gate_Rochester_large_opt1_10_1.5_no.7.qasm': continue #the solution of this circuit is corrupted
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
        
    