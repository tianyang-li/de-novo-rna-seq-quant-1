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
from libcpp.string cimport string

from misc_0 cimport ulong

cdef extern from "_graph_seq_0.h" namespace "_graph_seq_0":
    cdef cppclass SeqLoc:
        SeqLoc(ulong, ulong, ulong) except +
        
        ulong node_id
        ulong start
        ulong end
    
    ctypedef vector[SeqLoc] ReadNodeLoc
    
    cdef cppclass Node:
        Node(ulong, string) except +
        
        ulong node_id
        string seq_len
        vector[ulong] edges
    
    cdef cppclass PyGraph:
        PyGraph() except +
        
        ulong graph_id
        vector[Node] nodes
    
    cdef cppclass Fasta:
        Fasta() except +
        
        string info
        string seq 
    
    cdef cppclass ReadGraphLoc[T]:
        ReadGraphLoc(ulong) except +
            
        ulong graph_id
        vector[T] locs
    
    cdef cppclass ReadInGraph[T]:
        ReadInGraph() except +
        
        ulong read_id
        vector[ReadGraphLoc[T]] graph_locs


