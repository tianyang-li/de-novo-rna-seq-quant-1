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

from misc_0 cimport uint

cdef extern from "_graph_seq_0.h" namespace "_graph_seq_0":
    cdef cppclass SeqLoc:
        SeqLoc(uint, uint, uint) except +
        
        uint node_id
        uint start
        uint end
    
    ctypedef vector[SeqLoc] PyReadNodeLoc
    
    cdef cppclass PyNode:
        PyNode(uint, uint) except +
        
        uint node_id
        uint seq_len
        vector[uint] edges
    
    cdef cppclass PyGraph:
        PyGraph() except +
        
        uint graph_id
        vector[PyNode] nodes
    
    cdef cppclass Isoform:
        Isoform() except +
            
        uint graph_id


