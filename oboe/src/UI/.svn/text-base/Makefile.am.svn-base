# Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
# 
# OBOE is published under the Common Public License.
#
# Authors :
# Nidhi Sawhney <nsawhney@yahoo.com>
# The OBOE team
#
#

include_HEADERS = QpGenerator.h
#check_PROGRAMS = Test
#TESTS = Test

lib_LIBRARIES = libaccpm.a

libaccpm_a_SOURCES = QpGenerator.C QpGenerator.h

#Test_SOURCES = OracleTest.C
#Test_LDADD = libaccpmoracle.a  $(top_builddir)/src/AccpmLA/libaccpmla.a \
#	$(LAPACKCPP_DIR)/lib/liblapackpp.a -lblas -llapack
#Test_LDFLAGS = -L./ -L$(top_builddir)/src/AccpmLA -L$(LAPACKCPP_DIR)/lib/ 

INCLUDES = -I$(top_srcdir)/src/Utilities -I$(top_srcdir)/src/AccpmLA -I$(top_srcdir)/src/AccpmCore \
	-I$(top_srcdir)/src/ProblemInput -I$(top_srcdir)/src/Oracle \
	-I$(LAPACKCPP_DIR)/include/
