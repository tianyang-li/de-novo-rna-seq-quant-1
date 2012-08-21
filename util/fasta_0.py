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
some common FASTA related stuff
"""


class FastaSeq(object):
    __slots__ = ['seq', 'info']
    
    def __init__(self, info, seq):
        self.seq = seq
        self.info = info
    
    def __repr__(self):
        return self.info
    
    def __str__(self):
        return ">%s\n%s\n" % (self.info, self.seq)




