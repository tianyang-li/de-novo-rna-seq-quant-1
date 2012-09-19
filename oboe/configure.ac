# Copyright (C) 2005 Nidhi Sawhney <nsawhney@gmail.com>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([OBOE],
        [1.0],
        [Nidhi Sawhney nsawhney@gmail.com],
        [OBOE])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

AM_INIT_AUTOMAKE([dist-bzip2])
AC_PROG_INSTALL
AC_AIX
AC_ISC_POSIX
AC_MINIX
AC_STDC_HEADERS
dnl get the flags
dnl CFLAGS="${CFLAGS=}"
CXXFLAGS="${CXXFLAGS=}"

m4_include([m4/coin.m4])
m4_include([m4/oboe_usr_libs.m4])


dnl this macro is used to get the arguments supplied
dnl to the configure script (./configure --enable-debug)
dnl Check if we have enable debug support.
AC_MSG_CHECKING(whether to enable debugging)
default="no"
DEBUG=
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                       [default=$default]],, enable_debug=$default)
if test "x$enable_debug" = "xyes"; then
         CXXFLAGS="$CXXFLAGS -g -DDEBUG -Wall"
	 DEBUG=-g;
AC_MSG_RESULT(yes)
else
AC_MSG_RESULT(no)
fi

AC_MSG_CHECKING(whether to enable profiling)
default="no"
AC_ARG_ENABLE(debug, [  --enable-profile=[no/yes] turn on profiling
                       [default=$default]],, enable_profile=$default)
if test "x$enable_profile" = "xyes"; then
         CXXFLAGS="$CXXFLAGS -g -DDEBUG -pg"
	 DEBUG=-pg;
AC_MSG_RESULT(yes)
else
AC_MSG_RESULT(no)
fi

dnl And optimize, unless told otherwise
if test -z "$CXXFLAGS"; then CXXFLAGS="-O3 -Wall"; fi

AC_MSG_CHECKING(whether to enable serialization)
default="no"
AC_ARG_ENABLE(debug, [  --enable-serialization=[no/yes] turn on serialization support 
                       [default=$default]],, enable_serialization=$default)
if test "x$enable_serialization" = "xyes"; then
         CXXFLAGS="$CXXFLAGS -DSERIALIZATION"
 	 serialization=yes;
AC_MSG_RESULT(yes)
else
 	 serialization=no;
AC_MSG_RESULT(no)
fi

AC_COIN_HAS_USER_LIBRARY([Osi],[OSI],[OsiSolverInterface.hpp])
AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h])

dnl AC_OBOE_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],[glp_lpx_simplex])

dnl AC_PROG_CC
AC_PROG_CXX

dnl export these variable (so Makefile substitutions
dnl can be made. 
dnl AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(DEBUG)


dnl Libtool
dnl ---------------------------------------------------------------------- 
AM_PROG_LIBTOOL

AC_LANG_CPLUSPLUS
AC_PROG_RANLIB
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([for system type])
case "$host" in
	*-mingw32*)
	AC_DEFINE(WIN32, 1, [Define to 1 using WINDOWS])
	ARCH="WINDOWS"
	;;
	*-linux*)
	AC_DEFINE(LINUX, 1, [Define to 1 using LINUX])
	ARCH="LINUX"
	;;
esac
AC_MSG_RESULT("$ARCH")

dnl Look for blas - abandon all hope yea who enter here
blas_ok=yes
if test -z "$BLAS"; then
  blas_ok=no
  dnl The BLAS environment variable was not set
  dnl Look for blas in some typical places
  case "$host" in
  *-mingw32*)	
	BLAS=-lblas32
  ;;	
  sparc-sun-* )
    dnl AC_SEARCH_LIBS(dcopy,sunperf blas,,BLAS="libblas.a")
    AC_CHECK_LIB(sunperf,dcopy,BLAS=-lsunperf)
    echo "Checking for dcopy in -xlic_lib=sunperf (will look like I am"
    echo "testing -lm.)"
    if test -z "$BLAS"; then
      oldLIBS="$LIBS"
      LIBS="$LIBS -xlic_lib=sunperf"
      AC_CHECK_LIB(m,dcopy,BLAS=-xlic_lib=sunperf)
      LIBS="$oldLIBS"
    fi
    if test -z "$BLAS"; then
      AC_CHECK_LIB(blas,dcopy,BLAS=-lblas)
    fi
    dnl end of the tests for sun
    ;;
  *)
     AC_CHECK_LIB(blas,dcopy,BLAS=-lblas)
     ;;
  esac
  if test -n "$BLAS"; then
  dnl     Got it!
    blas_ok=yes
  fi
fi
dnl End the blas environment variable was not set
if test "$blas_ok" = no; then
    echo "ERROR: Cannot find BLAS library"
    echo "Set the environment variable BLAS to the BLAS library you want to use"
    echo "and run configure again"
    exit
fi
AC_SUBST(BLAS)
echo "Using BLAS library $BLAS"
dnl Done looking for blas

dnl Look for lapack 
lapack_ok=yes
if test -z "$LAPACK"; then
  lapack_ok=no
  dnl The LAPACK environment variable was not set
  dnl Look for lapack in some typical places
  case "$host" in
  *-mingw32*)	
	LAPACK=-llapack32
  ;;	
  sparc-sun-* )
    dnl AC_SEARCH_LIBS(dcopy,sunperf lapack,,LAPACK="liblapack.a")
    AC_CHECK_LIB(sunperf,dcopy,LAPACK=-lsunperf)
    echo "Checking for dcopy in -xlic_lib=sunperf (will look like I am"
    echo "testing -lm.)"
    if test -z "$LAPACK"; then
      oldLIBS="$LIBS"
      LIBS="$LIBS -xlic_lib=sunperf"
      AC_CHECK_LIB(m,dcopy,LAPACK=-xlic_lib=sunperf)
      LIBS="$oldLIBS"
    fi
    if test -z "$LAPACK"; then
      AC_CHECK_LIB(lapack,dcopy,LAPACK=-llapack)
    fi
    dnl end of the tests for sun
    ;;
  *)
     AC_CHECK_LIB(lapack,dcopy,LAPACK=-llapack)
     ;;
  esac
  if test -n "$LAPACK"; then
  dnl     Got it!
    lapack_ok=yes
  fi
fi
dnl End the lapack environment variable was not set
if test "$lapack_ok" = no; then
    echo "ERROR: Cannot find LAPACK library"
    echo "Set the environment variable LAPACK to the LAPACK library you want to use"
    echo "and run configure again"
    exit
fi
AC_SUBST(LAPACK)
echo "Using LAPACK library $LAPACK"
dnl Done looking for lapack

OBOE_LDFLAGS="";
dnl Check for LAPACKCPP
if test -z "$OBOE_ROOT"; then
OBOE_ROOT=`pwd`;
fi
AC_SUBST(OBOE_ROOT)

if test -z "$LAPACKCPP_DIR"; then
	LAPACKCPP_DIR="$OBOE_ROOT/lapackpp-2.1.0"	
	if test ! -d  "$LAPACKCPP_DIR"; then
		echo "" 	
		echo "ERROR: LAPACKCPP_DIR not found"
		echo "Please set the variable appropriately"
		echo "and run configure again"
		echo ""
		exit
	fi 	
fi

if test -z "$LAPACKCPP_DIR/include"; then
        echo ""
	echo "ERROR: LAPACKCPP_DIR does not have include directory"
	echo "Please set up LAPACKCPP_DIR such that LAPACKCPP_DIR/include"
	echo "or LAPACKCPP_DIR/include/lapackpp has the required header files."
	echo ""
	exit

fi  	

if test -z "$LAPACKCPP_LIB"; then
case "$host" in
  *-mingw32*)	
	LAPACKCPP_LIB=-llapackpp32
  ;;
  *)
	if test -d  "$LAPACKCPP_DIR/lib"; then
		OBOE_LDFLAGS="-L$LAPACKCPP_DIR/lib/"
		LAPACKCPP_LIB="$LAPACKCPP_DIR/lib/liblapackpp.a"
 		if test ! -f  "$LAPACKCPP_LIB"; then
		echo "" 	
		echo "ERROR: LAPACKCPP_LIB not found"
		echo "Please set the variable appropriately"
		echo "and run configure again"
		echo ""
		exit
		fi
	fi
   ;;
esac
fi

ACCPM_DIR=$OBOE_ROOT
AC_SUBST(ACCPM_DIR)
AC_SUBST(LAPACKCPP_DIR)
echo "Using LAPACK++ directory  $LAPACKCPP_DIR"
AC_SUBST(LAPACKCPP_LIB)
echo "Using LAPACK++ library  $LAPACKCPP_LIB"

if test "$serialization" = yes; then
OBOE_LDADD="-lboost_serialization $LAPACKCPP_LIB $LAPACK $BLAS -lg2c"
else
OBOE_LDADD="$LAPACKCPP_LIB $LAPACK $BLAS -lg2c"
fi

if test "$oboe_has_glpk" = true; then
OBOE_LDFLAGS="$OBOE_LDFLAGS -L$GLPKLIBDIR"
OBOE_LDADD="$OBOE_LDADD $GLPKLIBDIR/$GLPKLIB"
OBOE_INCLUDES="$OBOE_INCLUDES -I$GLPKINCDIR"
fi
if test "$coin_has_osi" = true; then
OBOE_LDFLAGS="$OBOE_LDFLAGS $OSILIB"
OBOE_INCLUDES="$OBOE_INCLUDES -I$OSIINCDIR"
fi

if test "$coin_has_glpk" = true; then
OBOE_LDFLAGS="$OBOE_LDFLAGS $GLPKLIB"
OBOE_INCLUDES="$OBOE_INCLUDES -I$GLPKINCDIR"
fi
AC_SUBST(OBOE_LDFLAGS)
AC_SUBST(OBOE_LDADD)
OBOE_INCLUDES="$OBOE_INCLUDES -I$LAPACKCPP_DIR -I$LAPACKCPP_DIR/include/ -I$LAPACKCPP_DIR/include/lapackpp"
AC_SUBST(OBOE_INCLUDES)


AC_CONFIG_FILES([
   Makefile
   README
   doc/Makefile
   m4/Makefile
   src/Makefile
   src/Utilities/Makefile
   src/AccpmLA/Makefile
   src/ProblemInput/Makefile
   src/AccpmCore/Makefile
   src/Oracle/Makefile
   src/UI/Makefile
   src/Benchmark/Makefile
   src/Benchmark/LA/Makefile
   src/Benchmark/PMED/Makefile
   src/Benchmark/QP/Makefile
   src/Benchmark/TSP/Makefile
   src/Benchmark/UAX/Makefile
   src/Benchmark/CS/Makefile
   src/Benchmark/HACCPM/Makefile
   src/Benchmark/NNPOL/Makefile
   src/Benchmark/MCF/Makefile
   src/Benchmark/APS/Makefile
])

AC_OUTPUT
