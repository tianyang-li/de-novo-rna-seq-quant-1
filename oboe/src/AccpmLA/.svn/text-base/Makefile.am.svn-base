# Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
# 
# OBOE is published under the Common Public License.
#
# Authors :
# Nidhi Sawhney <nsawhney@yahoo.com>
# The OBOE team
#
# 

include_HEADERS = AccpmVector.h AccpmGenMatrix.h AccpmDynMatrix.h AccpmBlasInterface.h

#check_PROGRAMS = Test TestSolve
check_PROGRAMS = Test

#Uncomment the following to run test with make check
#TESTS = $(check_PROGRAMS)

lib_LIBRARIES = libaccpmla.a

libaccpmla_a_SOURCES = AccpmVector.C AccpmVector.h AccpmGenMatrix.C AccpmGenMatrix.h  \
	AccpmDynMatrix.C  AccpmDynMatrix.h \
	AccpmLASolve.h AccpmLASolve.C AccpmBlasInterface.h AccpmBlasInterface.C

Test_SOURCES = MatrixTest.C 
Test_LDADD = libaccpmla.a $(OBOE_LDADD) 
Test_LDFLAGS = -L./ $(OBOE_LDFLAGS) 

#TestSolve_SOURCES = LinSolveTest.C
#TestSolve_LDADD = libaccpmla.a $(OBOE_LDADD) 
#TestSolve_LDFLAGS = -L./ $(OBOE_LDFLAGS) 

INCLUDES = -I$(top_srcdir)/src/Utilities \
	$(OBOE_INCLUDES)
