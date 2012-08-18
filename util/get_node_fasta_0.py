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

import networkx as nx


def get_node_fasta(dot_file, out_fasta):
    g = nx.read_dot(dot_file)
    with open(out_fasta, 'w') as fout:
        
        
    

def main():
    dot_file = None
    out_fasta = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'd:f:')
    except getopt.GetoptError as err:
        print >> sys.stderr, "missing"
        sys.exit(1)
    if (not dot_file
        or not out_fasta):
        print >> sys.stderr, "missing"
        sys.exit(1)

    get_node_fasta(dot_file, out_fasta)



