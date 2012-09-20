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

INCLUDES = -Igsl-1.15 -Isrc -Iboost_1_51_0 

LAPACK = /usr/lib/liblapack.a
BLAS = /usr/lib/libblas.a
LAPACKCPP_DIR = $(CURDIR)/lapackpp-2.5.4
LAPACKCPP_LIB = $(CURDIR)/build/lapackcpp/myinstall/lib/liblapackpp.a


all: util/single_1.so

util/single_1.so: lib/libgsl.a lib/oboe.a \
		setup.py lib/quant.a src/single_1.pyx  
	python setup.py build_ext -i 
	mv single_1.so util/
	
lib/oboe.a:
	tar xzvf oboe.tar.gz
	tar xzvf lapackpp-2.5.4.tar.gz
	mkdir -p build/lapackcpp/myinstall
	mkdir -p build/oboe/myinstall
	cd build/lapackcpp && \
		$(CURDIR)/lapackpp-2.5.4/configure \
		CFLAGS="-fPIC" CPPFLAGS="-fPIC" CXXFLAGS="-fPIC" \
		--prefix=$(CURDIR)/build/lapackcpp/myinstall && \
		make && make install 
	cd build/oboe && \
		$(CURDIR)/oboe/configure \
		CFLAGS="-fPIC" CXXFLAGS="-fPIC" CPPFLAGS="-fPIC" \
		FFLAGS="-fPIC" BLAS=$(BLAS) LAPACK=$(LAPACK) \
		LAPACKCPP_DIR=$(LAPACKCPP_DIR) LAPACKCPP_LIB=$(LAPACKCPP_LIB) \
		--prefix=$(CURDIR)/build/oboe/myinstall && \
		make && make install 
	mkdir -p build/oboe/objs && cd build/oboe/objs && \
		ar x ../myinstall/lib/libaccpm.a && \
		ar x ../myinstall/lib/libaccpmcore.a && \
		ar x ../myinstall/lib/libaccpmla.a && \
		ar x ../myinstall/lib/libaccpmoracle.a && \
		ar x ../myinstall/lib/libaccpmparam.a && \
		ar rcs $(CURDIR)/lib/oboe.a *.o

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
	#rm -rfv build/oboe #TODO make this into a command 
	#rm -rfv build/lapackcpp #TODO make this into a command 
	#rm -rfv gsl-1.15 #TODO make this into a command 
	#rm -rfv oboe #TODO make this into a command 
	#rm -rfv lapackpp-2.5.4 #TODO make this into a command 

.PHONY: all clean 


