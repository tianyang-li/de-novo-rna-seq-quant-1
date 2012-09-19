# Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
# 
# OBOE is published under the Common Public License.
#
# Authors :
# Nidhi Sawhney <nsawhney@yahoo.com>
# The OBOE team
#
#

include_HEADERS = Parameters.h

check_PROGRAMS = Test

#Uncomment the following to run test with make check
#TESTS = $(check_PROGRAMS)

lib_LIBRARIES = libaccpmparam.a

libaccpmparam_a_SOURCES = Parameters.C

Test_SOURCES = ParameterTest.C
Test_LDADD = libaccpmparam.a \
	$(top_builddir)/src/AccpmLA/libaccpmla.a \
 	$(OBOE_LDADD)	
Test_LDFLAGS = -L./ -L$(top_builddir)/src/AccpmLA $(OBOE_LDFLAGS) 
INCLUDES = -I$(top_srcdir)/src/AccpmLA/ -I$(top_srcdir)/src/Utilities $(OBOE_INCLUDES)

EXTRA_DIST = param.oboe
