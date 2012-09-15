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

from libc.stdio cimport printf

from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref, preincrement as inc

from graph_seq_0 cimport SeqLoc, ReadNodeLoc, PyGraph, Node, Fasta, ReadInGraph, ReadGraphLoc
from fasta_0 import FastaSeq
from misc_0 cimport ulong

cdef extern from "_single_1.h" namespace "_single_1":   
    ctypedef ReadNodeLoc SingleNodeLoc
    
    cdef void _get_isoforms(vector[PyGraph] * py_graphs,
                            vector[ReadInGraph[SingleNodeLoc]] * py_reads,
                            vector[Fasta] * isoforms,
                            ulong max_run)


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
    
    cdef ulong graph_id = 0
    
    cdef ulong node_id
    
    cdef Node * py_node 
    
    cdef vector[PyGraph].iterator py_graph_iter = py_graphs.begin()
    
    cdef vector[Node].iterator py_node_iter
    
    for graph_name, graph in graph_dict.iteritems():
        id_map = IdMap(graph_id)
        
        node_id = 0
        
        for node in graph.nodes():
            id_map.nodes_map[node] = node_id
            
            py_node = new Node(node_id,
                                 < string >< char *> graph.node[node]['label'])
            
            deref(py_graph_iter).nodes.push_back(deref(py_node))
            
            del py_node
            
            inc(py_node_iter)
            inc(node_id)
        
        py_node_iter = deref(py_graph_iter).nodes.begin() 
        
        for node in graph.nodes():
            for edge in graph[node]:
                deref(py_node_iter).edges.push_back(id_map.nodes_map[edge])
                
            inc(py_node_iter)
        
        id_maps[graph_name] = id_map 
        
        inc(graph_id)
        inc(py_graph_iter)


cdef void get_read_in_graph_from_py(read_in_graph, id_maps,
                                    vector[ReadInGraph[SingleNodeLoc]] * py_reads):        
    cdef ulong read_id = 0
    cdef ulong graph_id, node_id
    
    cdef vector[ReadInGraph[SingleNodeLoc]].iterator py_read_iter = py_reads.begin()
    
    cdef ReadGraphLoc[SingleNodeLoc] * py_graph_loc
    
    cdef ReadNodeLoc * py_node_loc
    
    cdef SeqLoc * seq_loc
    
    for read_info in read_in_graph.itervalues():
        deref(py_read_iter).read_id = read_id
        
        for graph_name, graph_aligns in read_info.iteritems():
            id_map = id_maps[graph_name]
            
            py_graph_loc = new ReadGraphLoc[SingleNodeLoc](id_map.graph_id)
            
            for graph_align in graph_aligns:
                py_node_loc = new ReadNodeLoc()
                
                for align in graph_align:
                    seq_loc = new SeqLoc(id_map.nodes_map[align.node],
                                         align.start, align.end)
                    
                    py_node_loc.push_back(deref(seq_loc))
                    
                    del seq_loc
                    
                py_graph_loc.locs.push_back(deref(py_node_loc))
                
                del py_node_loc 
            
            deref(py_read_iter).graph_locs.push_back(deref(py_graph_loc))
            
            del py_graph_loc
        
        inc(py_read_iter)
        inc(read_id)


cdef void convert_isoforms(isoforms, vector[Fasta] * _isoforms):
    """
    isoforms:
        a list of FastaSeq
    """
    
    cdef vector[Fasta].iterator isof_iter = _isoforms.begin()
    
    while isof_iter != _isoforms.end():
        isoforms.append(FastaSeq(deref(isof_iter).info.c_str(),
                                 deref(isof_iter).seq.c_str()))
        
        inc(isof_iter)
        

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
    
    cdef vector[ReadInGraph[SingleNodeLoc]] * py_reads = new vector[ReadInGraph[SingleNodeLoc]](len(read_in_graph))
    
    get_read_in_graph_from_py(read_in_graph, id_maps, py_reads)
    
    cdef vector[Fasta] * _isoforms = new vector[Fasta]()

    _get_isoforms(py_graphs, py_reads, _isoforms, max_run)

    isoforms = []
    
    convert_isoforms(isoforms, _isoforms)

    del py_graphs
    del py_reads
    del _isoforms

    return isoforms


