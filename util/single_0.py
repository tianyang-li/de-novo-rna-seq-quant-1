#!/usr/bin/env python

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


import getopt
import sys

from Bio import SeqIO
import networkx as nx

from blat_0 import read_psl


class ReadInGraph(object):
    __slots__ = ['graph', 'node', 'loc']
    
    def __init__(self, graph, node, loc):
        # graph name
        self.graph = graph
        
        # node name
        self.node = node
        
        # start location of read on node
        self.loc = loc


def read_all_in_graph(read_graph):
    """
    read_graph 
        a dictionary containing list of where 
        the read is in the graph
    
    locate all the reads in the graph that lie
    completely within a node
    """


def main():
    dot_file = None
    
    # read -> contig 
    rcont_psl_file = None  
    
    # read -> each node's sequence
    rcomp_psl_file = None
    
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'd:', ['rcont='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    
    for opt, arg in opts:
        if opt == '-d':
            dot_file = arg
        if opt == '--rc':
            rcont_psl_file = arg
    
    if (not dot_file
        or not rcont_psl_file):
        print >> sys.stderr, "missing"
        sys.exit(1)
        

if __name__ == '__main__':
    main()        
        



