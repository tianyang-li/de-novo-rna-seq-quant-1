# Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
#  
# OBOE is published under the Common Public License.
#
# Authors :
# Nidhi Sawhney <nsawhney@yahoo.com>
# The OBOE team
# 
#

check_PROGRAMS = Test

#Uncomment the following to run test with make check
#TESTS = $(check_PROGRAMS)

lib_LIBRARIES = libaccpmcore.a

libaccpmcore_a_SOURCES = Manager.C Manager.h Method.C Method.h LocSet.C LocSet.h Solution.h Solution.C

Test_SOURCES = ManagerTest.C
Test_LDADD = libaccpmcore.a \
	$(top_builddir)/src/ProblemInput/libaccpmparam.a \
	$(top_builddir)/src/AccpmLA/libaccpmla.a \
 	$(OBOE_LDADD)	

Test_LDFLAGS = -L./ -L$(top_builddir)/src/ProblemInput \
	-L$(top_builddir)/src/AccpmLA \
	$(OBOE_LDFLAGS)	

INCLUDES = -I$(top_srcdir)/src/ProblemInput \
	-I$(top_srcdir)/src/AccpmLA/ \
	-I$(top_srcdir)/src/Oracle \
	-I$(top_srcdir)/src/Utilities \
	$(OBOE_INCLUDES)

EXTRA_DIST = param.txt
