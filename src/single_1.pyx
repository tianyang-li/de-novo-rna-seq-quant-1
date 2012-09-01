#  Copyright (C) 2012 Tianyang Li
#  tmy1018@gmail.com
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License

from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as inc

from graph_seq_0 cimport SeqLoc
from fasta_0 import FastaSeq
from misc_0 cimport uint

cdef extern from "_single_1.h" namespace "_single_1":
    cdef cppclass PyNode:
        PyNode(uint) except +
        
        uint node_id
        vector[uint] edges
    
    cdef cppclass PyGraph:
        PyGraph() except +
        
        uint graph_id
        vector[PyNode] nodes


class IdMap(object):
    __slots__ = ['graph_id', 'nodes_map']
    
    def __init__(self, graph_id):
        self.graph_id = graph_id
        self.nodes_map = {}
    
    def __str__(self):
        return "graph:%s,nodes_map:%s" % (self.graph_id, str(self.nodes_map))


cdef void get_graph_from_py(graph_dict, id_maps, vector[PyGraph] * py_graphs):
    """
    id_maps: 
        a dictionary where
    
            id_maps[graph] 
        
        is a IdMap that says where the graph is in py_graphs
    """
    
    cdef uint graph_id = 0
    cdef uint node_id
    cdef PyNode * py_node 
    
    cdef vector[PyGraph].iterator py_graph_iter = py_graphs.begin()
    cdef vector[PyNode].iterator py_node_iter
    
    for graph_name, graph in graph_dict.iteritems():
        id_map = IdMap(graph_id)
        
        node_id = 0
        
        for node in graph.nodes():
            id_map.nodes_map[node] = node_id
            py_node = new PyNode(node_id)
            deref(py_graph_iter).nodes.push_back(deref(py_node))
            del py_node
            #TODO: no dynamic mem allocation? faster?
            
            inc(py_node_iter)
            inc(node_id)
        
        py_node_iter = deref(py_graph_iter).nodes.begin() 
        
        for node in graph.nodes():
            for edge in graph[node]:
                deref(py_node_iter).edges.push_back(id_map.nodes_map[edge])
        
        id_maps[graph_name] = id_map 
        
        inc(graph_id)
        inc(py_graph_iter)
        

def get_isoforms(read_in_graph, graph_dict, max_run=1000000):
    """
    read_in_graph: 
        same as in single_0.py
    
    graph_dict: 
        same as in single_0.py
    
    max_run: 
        max number of runs in MCMC
    
    return the isoforms as a list of FastaSeq
    """
    
    cdef vector[PyGraph] * py_graphs = new vector[PyGraph](len(graph_dict))
    
    id_maps = {}
    
    get_graph_from_py(graph_dict, id_maps, py_graphs)
    
    isoforms = []
    
    del py_graphs
    
    return isoforms

