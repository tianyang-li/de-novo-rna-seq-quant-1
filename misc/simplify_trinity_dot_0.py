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
import getopt

import networkx as nx

def simplify_trinity_dot(in_file, out_file):
    g = nx.read_dot(in_file)

def main():
    in_file = None
    out_file = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'i:o:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-i':
            in_file = arg
        if opt == '-o':
            out_file = arg
    if (not in_file
        or not out_file):
        print >> sys.stderr, "missing"
        sys.exit(1)
        
    simplify_trinity_dot(in_file, out_file)

