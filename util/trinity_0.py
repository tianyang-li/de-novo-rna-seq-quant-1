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

from collections import defaultdict

from Bio import SeqIO

from fasta_0 import FastaSeq


def get_contig_dict(trinity_out_file):
    """
    return a defaultdict(dict)
    where
        
        contig_dict[graph][contig]
    
    is __contig__ from __graph__
        
    """
    
    contig_dict = defaultdict(dict)
    
    for rec in SeqIO.parse(trinity_out_file, 'fasta'):
        rec_id = rec.id
        contig_dict[rec_id.split("_")[0]][rec_id] = FastaSeq(rec_id,
                                                             str(rec.seq))
    
    return contig_dict


