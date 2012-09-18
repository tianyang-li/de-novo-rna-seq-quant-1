CC = g++

CFLAGS = -Wall -fPIC -Wconversion -Wextra -ggdb -DDEBUG \
-Wall -W -Wmissing-prototypes \
-Wstrict-prototypes \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common \
#-Wshadow \
-ansi -pedantic \
-Wpointer-arith -Wcast-qual -Wcast-align 
 
AR = ar

INCLUDES = -Igsl-1.15 -Isrc -Iboost_1_51_0 -Inlopt-2.3/myinstall/include

all: util/single_1.so

util/single_1.so: lib/libgsl.a lib/libnlopt_cxx.a \
		setup.py lib/quant.a src/single_1.pyx  
	python setup.py build_ext -i 
	mv single_1.so util/

src/single_1.pyx: src/graph_seq_0.pxd

src/graph_seq_0.pxd: src/misc_0.pxd src/_graph_seq_0.h src/_graph_seq_0.cc

src/misc_0.pxd: src/_misc_0.cc src/_misc_0.h

lib/libnlopt_cxx.a:
	tar xzvf nlopt-2.3.tar.gz
	mkdir nlopt-2.3/myinstall
	cd nlopt-2.3 && ./configure --with-cxx CXXFLAGS="-fPIC" \
		CFLAGS="-fPIC" CPPFLAGS="-fPIC" \
		--prefix=$(CURDIR)/nlopt-2.3/myinstall && \
		make && make install
	cp nlopt-2.3/myinstall/lib/libnlopt_cxx.a lib/libnlopt_cxx.a

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
	#rm -rfv nlopt-2.3 #TODO make this into a command 
	#rm -rfv gsl-1.15 #TODO make this into a command 

learn:
	cd learn && $(MAKE)

.PHONY: all clean learn


