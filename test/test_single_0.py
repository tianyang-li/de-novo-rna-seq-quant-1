#!/usr/bin/env python

import getopt
import sys


from util.single_0 import get_contig_dict, get_rcc_psl


def main():
    rcont_psl_file = None
    contig_file = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'c:', ['rcont='])
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-c':
            contig_file = arg
        if opt == '--rcont':
            rcont_psl_file = arg
    if (not rcont_psl_file
        or not contig_file):
        print >> sys.stderr, "missing"
        sys.exit(1)
    
    rcc_psl = get_rcc_psl(rcont_psl_file)
    
    contig_dict = get_contig_dict(contig_file)
    
    for comp_name, comp in rcc_psl.iteritems():
        for cont_name, psls in comp.iteritems():
            for psl in psls:
                node_start = contig_dict[comp_name][cont_name].find_start(psl.tStart)
                if (psl.tStart < node_start[1] 
                    or psl.tStart >= node_start[2]):
                    print >> sys.stderr, "error"
                node_end = contig_dict[comp_name][cont_name].find_end(psl.tEnd)
                if (psl.tEnd <= node_end[1]
                    or psl.tEnd > node_end[2]):
                    print >> sys.stderr, "error"
    

if __name__ == '__main__':
    main()


