CC = g++
CFLAGS = -Wall
AR = ar
INCLUDES = -Igsl-1.15 -Isrc -Iboost_1_51_0

all: util/single_1.so

util/single_1.so: src/single_1.pyx lib/quant.a setup.py lib/libgsl.a 
	python setup.py build_ext -i
	mv single_1.so util/

src/single_1.pyx: src/graph_seq_0.pxd

src/graph_seq_0.pxd: src/misc_0.pxd src/_graph_seq_0.h src/_graph_seq_0.cc

src/misc_0.pxd: src/_misc_0.cc src/_misc_0.h

lib/libgsl.a:
	tar xzvf gsl-1.15.tar.gz
	cd gsl-1.15 && ./configure && make
	cp gsl-1.15/.libs/libgsl.a lib/

lib/quant.a: build/_single_1.o build/_graph_seq_0.o build/_misc_0.o
	$(AR) rcs lib/quant.a build/_single_1.o build/_graph_seq_0.o \
		build/_misc_0.o

build/_misc_0.o: src/_misc_0.cc src/_misc_0.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/_misc_0.cc -o build/_misc_0.o

build/_single_1.o: src/_single_1.cc src/_single_1.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/_single_1.cc -o build/_single_1.o

build/_graph_seq_0.o: src/_graph_seq_0.cc src/_graph_seq_0.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/_graph_seq_0.cc -o build/_graph_seq_0.o
	
clean:
	rm -rfv build/*
	rm -rfv util/single_1.so
	rm -rfv lib/*
	rm -rfv src/*.cpp
	rm -rfv gsl-1.15/

.PHONY: all clean


