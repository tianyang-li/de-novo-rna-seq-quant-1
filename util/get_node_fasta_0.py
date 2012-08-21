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

import sys

import networkx as nx


def get_node_fasta(dot_files):
    for dot_file in dot_files:
        
        g = nx.read_dot(dot_file)
        
        name = g.graph['name']
        
        for n, n_info in g.node.iteritems():
            print ">%s:%s\n%s\n" % (name, n,
                                    n_info['label'].split("(")[0]),


def main():
    if not sys.argv[1:]:
        print >> sys.stderr, "missing"

    get_node_fasta(sys.argv[1:])


if __name__ == '__main__':
    main()
