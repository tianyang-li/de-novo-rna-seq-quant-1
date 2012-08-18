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


def get_cr_dict(rc_psl_file):
    """
    rc_psl_file
        sorted by tName then by tStart
    """
    
    cur_tName = None
    cur_cr_list = None

    cr_dict = {}
    
    for psl in read_psl(rc_psl_file):
        if psl.tName != cur_tName:
            cr_dict[cur_tName] = cur_cr_list
            cur_tName = psl.tName
            cur_cr_list = []
        cur_cr_list.append(psl)
    
    del cr_dict[None]

    return cr_dict


def main():
    dot_file = None
    
    # sorted by tName then by tStart
    rc_psl_file = None  
    
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'd:', ['rc='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    
    for opt, arg in opts:
        if opt == '-d':
            dot_file = arg
        if opt == '--rc':
            rc_psl_file = arg
    
    if (not dot_file
        or not rc_psl_file):
        print >> sys.stderr, "missing"
        sys.exit(1)
    
    cr_dict = get_cr_dict(rc_psl_file)
        

if __name__ == '__main__':
    main()        
        



