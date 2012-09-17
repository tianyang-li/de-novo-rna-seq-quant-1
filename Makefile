CC = g++

CFLAGS = -Wall -fPIC -Wconversion -Wextra -ggdb -DDEBUG \
-DHAVE_STD -DHAVE_NAMESPACES \
-Wall -W -Wmissing-prototypes \
-Wstrict-prototypes \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common \
#-Wshadow \
-ansi -pedantic \
-Wpointer-arith -Wcast-qual -Wcast-align 
 
AR = ar

INCLUDES = -Igsl-1.15 -Isrc -Iboost_1_51_0 -Ioptpp-2.4/include \
-Ioptpp-2.4/newmat11

all: util/single_1.so

util/single_1.so: lib/libnewmat.a lib/libopt.a lib/libgsl.a setup.py \
	lib/quant.a src/single_1.pyx  

	python setup.py build_ext -i 
	mv single_1.so util/

lib/libopt.a:
	tar xzvf optpp-2.4.tar.gz
	
	cd optpp-2.4 && ./configure CFLAGS="-fPIC" CPPFLAGS="-fPIC" \
		FFLAGS="-fPIC" && $(MAKE)
	
	cp optpp-2.4/lib/.libs/*.a lib

lib/libnewmat.a:
	tar xzvf optpp-2.4.tar.gz
	
	cd optpp-2.4 && ./configure CFLAGS="-fPIC" CPPFLAGS="-fPIC" \
		FFLAGS="-fPIC" && $(MAKE)
	
	cp optpp-2.4/lib/.libs/*.a lib

src/single_1.pyx: src/graph_seq_0.pxd

src/graph_seq_0.pxd: src/misc_0.pxd src/_graph_seq_0.h src/_graph_seq_0.cc

src/misc_0.pxd: src/_misc_0.cc src/_misc_0.h

lib/libgsl.a:
	tar xzvf gsl-1.15.tar.gz
	cd gsl-1.15 && ./configure && $(MAKE)
	cp gsl-1.15/.libs/libgsl.a lib/

lib/quant.a: build/_single_1.o build/_graph_seq_0.o build/_misc_0.o \
				build/_mcmc_0.o
	$(AR) rcs lib/quant.a build/_single_1.o build/_graph_seq_0.o \
		build/_misc_0.o build/_mcmc_0.o

build/_misc_0.o: src/_misc_0.cc src/_misc_0.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/_misc_0.cc -o build/_misc_0.o

build/_single_1.o: src/_single_1.cc src/_single_1.h src/_graph_seq_0.h src/_mcmc_0.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/_single_1.cc -o build/_single_1.o

build/_graph_seq_0.o: src/_graph_seq_0.cc src/_graph_seq_0.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/_graph_seq_0.cc -o build/_graph_seq_0.o

build/_mcmc_0.o: src/_mcmc_0.h src/_mcmc_0.cc src/_graph_seq_0.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/_mcmc_0.cc -o build/_mcmc_0.o
	
clean:
	rm -rfv build/*
	rm -rfv util/single_1.so
	rm -rfv lib/*
	rm -rfv src/*.cpp
	#rm -rfv gsl-1.15/ #TODO make this into a command 

learn:
	cd learn && $(MAKE)

.PHONY: all clean learn


