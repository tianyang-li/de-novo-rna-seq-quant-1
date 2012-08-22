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
from collections import defaultdict

from blat_0 import read_psl
from dot_0 import get_splice_graph
from trinity_0 import get_contig_dict


class NodeSeq(object):
    """
    specifies where a sequence is on
    in a graph's node
    
    it is assumed that 
    the graph this node is in
    is already known
    """
    
    __slots__ = ['node', 'start', 'end']
    
    def __init__(self, node, start, end):
        """
        start, end - Python index 
        """
        
        self.node = node
        self.start = start
        self.end = end


class ReadInGraph(object):
    __slots__ = ['graph', 'node_seqs']
    
    def __init__(self, graph, node, loc):
        # graph name
        self.graph = graph
        
        # a list of NodeSeq saying where the
        # read is in the graph
        self.node_seqs = []
    

def get_cnr_psl(rcomp_psl_file):
    """
    return dictionary where
    
        cnr_psl[graph][node]
    
    is a list of read aligned to that node 
    
    c - comp
    n - node
    r - read
    """
    
    cnr_psl = defaultdict(lambda : defaultdict(list))
    
    for psl in read_psl(rcomp_psl_file):
        comp, node = psl.tName.split(":")
        cnr_psl[comp][node].append(psl)
    
    return cnr_psl


def get_ccr_psl(rcont_psl_file):
    """
    return a dictionary where
    
        ccr_psl[graph]
    
    is a list of read aligned to a contig 
    in that graph
    
    1st c - comp
    2nd c - contig
    r - read
    """
    
    ccr_psl = defaultdict(lambda : defaultdict(list))
    
    for psl in read_psl(rcont_psl_file):
        tName = psl.tName
        ccr_psl[tName.split("_")[0]][tName].append(psl)
    
    return ccr_psl


def read_all_in_node(psl):
    """
    max allowed gap length is 1 bp
    """
    
    if psl.qSize - psl.qEnd > 1:
        return False
    if (psl.qSize - psl.qEnd == 1 and
        psl.tEnd == psl.tSize):
        return False 
    
    if psl.qStart > 1:
        return False
    if (psl.qStart == 1 and
        psl.tStart == 0):
        return False
    
    for i in xrange(psl.blockCount - 1):
        if psl.qStarts[i] + psl.blockSizes[i] + 1 < psl.qStarts[i + 1]:
            return False
        if psl.tStarts[i] + psl.blockSizes[i] + 1 < psl.tStarts[i + 1]:
            return False
    
    #TODO: have I missied some cases?
    
    return True


def get_rcc_psl(rcont_psl_file):


def read_all_in_graph(read_in_graph, cnr_psl,
                      splice_graph):
    """
    read_in_graph 
        a defaultdict(lambda : defaultdict(list)) 
        containing list of where 
        the read is in the graph
            
            read_in_graph[read][graph]
        
        is a list of alignments of read to nodes in 
        that graph
    
    cnr_psl:
        a dictionary where 
            cnr_psl[graph_name][node_name]
        is a list of psl alignment to that node
    
    locate all the reads in the graph that lie
    completely within a node
    """
    
    comp = splice_graph.name
    node_psl = cnr_psl[comp]
    
    for node in splice_graph.node:
        for psl in node_psl[node]:
            if read_all_in_node(psl):
                read_in_graph[psl.qName][comp].append([NodeSeq(node, psl.tStart - psl.qStart, psl.tEnd + psl.tSize - psl.tEnd)])


def read_across_node(read_in_graph, ccr_psl, contig_dict):
    """
    read_in_graph 
        a defaultdict(lambda : defaultdict(list)) 
        containing list of where 
        the read is in the graph
            
            read_in_graph[read][graph]
        
        is a list of alignments of read to nodes in 
        that graph
    """
                    

def main():
    dot_file = None
    
    # read -> contig 
    rcont_psl_file = None  
    
    # read -> each node's sequence
    rcomp_psl_file = None

    contig_file = None
    
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'd:c:',
                                ['rcont=', 'rcomp='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    
    for opt, arg in opts:
        if opt == '-d':
            dot_file = arg
        if opt == '--rcont':
            rcont_psl_file = arg
        if opt == '--rcomp':
            rcomp_psl_file = arg
        if opt == '-c':
            contig_file = arg
    
    if (not dot_file
        or not contig_file
        or not rcont_psl_file
        or not rcomp_psl_file):
        print >> sys.stderr, "missing"
        sys.exit(1)
    
    read_in_graph = defaultdict(lambda : defaultdict(list))
    
    cnr_psl = get_cnr_psl(rcomp_psl_file)
    
    splice_graph = get_splice_graph(dot_file)
    
    read_all_in_graph(read_in_graph, cnr_psl, splice_graph)
    
    contig_dict = get_contig_dict(contig_file)
    
    print >> sys.stderr, "hello world"
        

if __name__ == '__main__':
    main()        
        



