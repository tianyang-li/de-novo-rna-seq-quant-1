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
matches - Number of bases that match that aren't repeats
misMatches - Number of bases that don't match
repMatches - Number of bases that match but are part of repeats
nCount - Number of 'N' bases
qNumInsert - Number of inserts in query
qBaseInsert - Number of bases inserted in query
tNumInsert - Number of inserts in target
tBaseInsert - Number of bases inserted in target
strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
qName - Query sequence name
qSize - Query sequence size
qStart - Alignment start position in query
qEnd - Alignment end position in query
tName - Target sequence name
tSize - Target sequence size
tStart - Alignment start position in target
tEnd - Alignment end position in target
blockCount - Number of blocks in the alignment (a block contains no gaps)
blockSizes - Comma-separated list of sizes of each block
qStarts - Comma-separated list of starting positions of each block in query
tStarts - Comma-separated list of starting positions of each block in target 
"""

import sys

from itertools import izip

class BlatEntry(object):
    __slots__ = ["matches", "misMatches",
                 "repMatches", "nCount",
                 "qNumInsert", "qBaseInsert",
                 "tNumInsert", "tBaseInsert",
                 "strand", "qName", "qSize",
                 "qStart", "qEnd", "tName",
                 "tSize", "tStart", "tEnd",
                 "blockCount", "blockSizes",
                 "qStarts", "tStarts"]
    
    def __init__(self, line):
        entries = line.strip().split("\t")
        self.matches = int(entries[0])
        self.misMatches = int(entries[1])
        self.repMatches = int(entries[2])
        self.nCount = int(entries[3])
        self.qNumInsert = int(entries[4])
        self.qBaseInsert = int(entries[5])
        self.tNumInsert = int(entries[6])
        self.tBaseInsert = int(entries[7])
        self.strand = entries[8]
        self.qName = entries[9]
        self.qSize = int(entries[10])
        self.qStart = int(entries[11])
        self.qEnd = int(entries[12])
        self.tName = entries[13]
        self.tSize = int(entries[14])
        self.tStart = int(entries[15])
        self.tEnd = int(entries[16])
        self.blockCount = int(entries[17])
        self.blockSizes = map(lambda s: int(s), entries[18].split(",")[:-1])
        self.qStarts = map(lambda s: int(s), entries[19].split(",")[:-1])
        self.tStarts = map(lambda s: int(s), entries[20].split(",")[:-1])
    
    def __str__(self):
        return (("%d\t" % self.matches)
                + ("%d\t" % self.misMatches)
                + ("%d\t" % self.repMatches)
                + ("%d\t" % self.nCount)
                + ("%d\t" % self.qNumInsert)
                + ("%d\t" % self.qBaseInsert)
                + ("%d\t" % self.tNumInsert)
                + ("%d\t" % self.tBaseInsert)
                + self.strand + "\t"
                + self.qName + "\t"
                + ("%d\t" % self.qSize)
                + ("%d\t" % self.qStart)
                + ("%d\t" % self.qEnd)
                + self.tName + "\t"
                + ("%d\t" % self.tSize)
                + ("%d\t" % self.tStart)
                + ("%d\t" % self.tEnd)
                + ("%d\t" % self.blockCount)
                + "".join(("%d," % x) for x in self.blockSizes) + "\t"
                + "".join(("%d," % x) for x in self.qStarts) + "\t"
                + "".join(("%d," % x) for x in self.tStarts))
    
    def get_tBlocks(self):
        return self._get_blocks(self.tStarts)
    
    def get_qBlocks(self):
        return self._get_blocks(self.qStarts)
        
    def _get_blocks(self, starts):
        """
        each block is a Pythonic range
        """
        blocks = []
        for start, block_size in izip(starts, self.blockSizes):
            blocks.append((start, start + block_size))
        

def read_psl(psl_file):
    """
    for
        SINGLE
    reads only
    """
    with open(psl_file, 'r') as fin:
        for line in fin:
            yield BlatEntry(line)


def main():
    if not sys.argv[1:]:
        print >> sys.stderr, "missing"
    for arg in sys.argv[1:]:
        for entry in read_psl(arg):
            print str(entry)

if __name__ == '__main__':
    main()



