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
    
    name = g.graph['name']
    
    with open(out_fasta, 'w') as fout:
        for n, n_info in g.node.iteritems():
            fout.write(">%s:%s\n%s\n" % 
                       (name, n,
                        n_info['label'].split("(")[0]))
    

def main():
    dot_file = None
    out_fasta = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'd:f:')
    except getopt.GetoptError as err:
        print >> sys.stderr, "missing"
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-d':
            dot_file = arg
        if opt == '-f':
            out_fasta = arg
    if (not dot_file
        or not out_fasta):
        print >> sys.stderr, "missing"
        sys.exit(1)

    get_node_fasta(dot_file, out_fasta)


if __name__ == '__main__':
    main()
