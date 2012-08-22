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

"""
handles Trinity's output
"""

from __future__ import division

from collections import defaultdict
import re

from Bio import SeqIO

from fasta_0 import FastaSeq


class TrinityContig(FastaSeq):
    __slots__ = ['nodes']
    
    def __init__(self, rec_id, seq, nodes):
        """
        nodes
            a list of
             
                (node, low, high)
            
            describing nodes in splice graph
        """
        
        super(TrinityContig, self).__init__(rec_id, seq)
        
        self.nodes = nodes
    
    def find_start(self, start):
        """
        find the node of a psl.tStart is in
        
        it is assumed that _start_ can always
        find a node to fit in
        """
        l = 0
        h = len(self.nodes) - 1
        x = int((l + h) / 2)
        while (start < self.nodes[x][1]
               or start >= self.nodes[x][2]):
            if start < self.nodes[x][1]:
                h = x - 1
            else:
                l = x + 1
            x = int((l + h) / 2)  
        return self.nodes[x]
    
    def find_end(self, end):
        """
        find the node a psl.tEnd is in
        
        it is assumed that _end_ can always find 
        a node to fit in
        """
        l = 0
        h = len(self.nodes) - 1
        x = int((l + h) / 2)
        while (end <= self.nodes[x][1]
               or end > self.nodes[x][2]):
            if end <= self.nodes[x][1]:
                h = x - 1
            else:
                l = x + 1
            x = int((l + h) / 2)
        return self.nodes[x]


path_re = re.compile(r'path=\[(.*)\]')
node_re = re.compile(r'(\w*):(\d*)-(\d*)')

def convert_node_re(n):
    return (n[0], int(n[1]), int(n[2]) + 1)

def get_contig_nodes(rec_decription):
    global path_re, node_re
    
    nodes = path_re.search(rec_decription).group(1).split(" ")
    
    return map(lambda n: convert_node_re(node_re.search(n).groups()), nodes)


def get_contig_dict(trinity_out_file):
    """
    return a defaultdict(dict)
    where
        
        contig_dict[graph][contig]
    
    is _contig_ from _graph_
        
    """
    
    contig_dict = defaultdict(dict)
    
    for rec in SeqIO.parse(trinity_out_file, 'fasta'):
        rec_id = rec.id
        contig_dict[rec_id.split("_")[0]][rec_id] = TrinityContig(rec_id, str(rec.seq), get_contig_nodes(rec.description))
    
    return contig_dict


