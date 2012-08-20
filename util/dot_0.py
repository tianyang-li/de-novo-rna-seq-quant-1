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
some functions that are required with the dot files
"""

import networkx as nx

def get_splice_graph(dot_file):
    g = nx.read_dot(dot_file)
    
    for n in g.node:
        g.node[n] = {'label': g.node[n]['label'].split("(")[0]}
    
    for n1, n2 in g.edges():
        g[n1][n2][0] = {}
        
    # a list of reads in this graph
    g.graph['reads'] = []
    
    return g


