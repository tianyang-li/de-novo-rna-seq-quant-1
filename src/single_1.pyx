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
        

def get_isoforms(read_in_graph, graph_dict):
    """
    read_in_graph: same as in single_0.py
    
    graph_dict: same as in single_0.py
    
    return the isoforms as a list of FastaSeq
    """
    
    cdef vector[PyGraph] *py_graphs = new vector[PyGraph](len(graph_dict))
    
    isoforms = []
    
    del py_graphs
    
    return isoforms

