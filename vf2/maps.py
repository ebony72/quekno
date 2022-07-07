#-*- coding:utf-8 -*-
# AUTHOR:   yaolili
# FILE:     map.py
# ROLE:     current map result
# CREATED:  2015-11-28 20:55:11
# MODIFIED: 2015-12-04 16:01:43
# ADDAPTED: 2019-06-25 for Quantum Circuit Transformation by Sanjiang Li (SL)
# all comments by SL started with '#sl'

import networkx as nx

class Map:   
    #sl result is a dict, i.e., a mapping from logical qubits to physical qubits
    def __init__(self, result):  
        self.__subMap = []
        self.__gMap = []    
        if type(result) is not dict:
            print("Class Map __init__() argument type error! dict expected!")
            exit()
        if result:             
            for key in result:     
                self.__subMap.append(key)
                self.__gMap.append(result[key])
                
    def subMap(self): #sl: used when checking isMeetRules 
        return self.__subMap
        
    def gMap(self):
        return self.__gMap
            
    #type = 0, g=subGraph; type = 1, g=graph
    def neighbor(self, g, type):

        if not (type == 1 or type == 0):
            print("Class Map neighbor() argument value error! type expected 0 or 1!")
            exit()
                    
        if type:
            curMap = self.__gMap
        else:
            curMap = self.__subMap

        neighbor_set = set()
        for x in curMap:
            for q in nx.all_neighbors(g,x):
                neighbor_set.add(q)
        '''Return unmapped neighbors''' #sl
        neighbor = list(neighbor_set - set(curMap))
        if not neighbor and not neighbor_set <= set(curMap):
            raise Exception('Neighbor cannot be empty!')
                    
        return neighbor  
