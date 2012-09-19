# Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
# 
# OBOE is published under the Common Public License.
#
# Authors :
# Nidhi Sawhney <nsawhney@yahoo.com>
# The OBOE team
#
#

check_PROGRAMS = oboeTSP 
#TESTS = oboeTSP 

oboeTSP_SOURCES = TSPTest.C
oboeTSP_LDADD = $(top_builddir)/src/UI/libaccpm.a \
		$(top_builddir)/src/AccpmCore/libaccpmcore.a \
		$(top_builddir)/src/AccpmLA/libaccpmla.a \
		$(top_builddir)/src/ProblemInput/libaccpmparam.a \
		$(top_builddir)/src/Oracle/libaccpmoracle.a \
		$(OBOE_LDADD)

oboeTSP_LDFLAGS = -L./ -L$(top_builddir)/src/UI \
	-L$(top_builddir)/src/AccpmCore \
	-L$(top_builddir)/src/AccpmLA \
	-L$(top_builddir)/src/ProblemInput \
	-L$(top_builddir)/src/Oracle \
	$(OBOE_LDFLAGS) 

INCLUDES = -I$(top_srcdir)/src/Utilities \
	 -I$(top_srcdir)/src/AccpmLA \
	-I$(top_srcdir)/src/AccpmCore \
	-I$(top_srcdir)/src/ProblemInput \
	-I$(top_srcdir)/src/Oracle \
	-I$(top_srcdir)/src/UI \
	$(OBOE_INCLUDES)

EXTRA_DIST = param.txt tsp29.coord
