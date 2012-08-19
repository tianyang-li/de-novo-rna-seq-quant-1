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


def get_node_fasta(dot_files, out_file):
    with open(out_file, 'w') as fout:
        for dot_file in dot_files:
            
            g = nx.read_dot(dot_file)
            
            name = g.graph['name']
            
            for n, n_info in g.node.iteritems():
                fout.write(">%s:%s\n%s\n" % 
                           (name, n,
                            n_info['label'].split("(")[0]))
    

def main():
    out_file = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], '-o')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    if not args or not out_file:
        print >> sys.stderr, "missing"
        sys.exit(1)

    get_node_fasta(args, out_file)


if __name__ == '__main__':
    main()
