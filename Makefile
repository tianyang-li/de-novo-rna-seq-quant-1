CC = g++
CFLAGS = -Wall
AR = ar

all: util/single_1.so

util/single_1.so: src/single_1.pyx lib/quant.a setup.py lib/libgsl.a
	python setup.py build_ext -i
	mv single_1.so util/

lib/libgsl.a:
	tar xzvf gsl-1.15.tar.gz
	cd gsl-1.15 && ./configure && make
	cp gsl-1.15/.libs/libgsl.a lib/

lib/quant.a: build/single_1.o build/graph_seq_0.o
	$(AR) rcs lib/quant.a build/single_1.o build/graph_seq_0.o

build/single_1.o: src/single_1.cc src/single_1.h
	$(CC) $(CFLAGS) -c src/single_1.cc -o build/single_1.o

build/graph_seq_0.o: src/graph_seq_0.cc src/graph_seq_0.h
	$(CC) $(CFLAGS) -c src/graph_seq_0.cc -o build/graph_seq_0.o
	
clean:
	rm -rfv build/*
	rm -rfv util/single_1.so
	rm -rfv lib/*
	rm -rfv src/single_1.cpp
	rm -rfv gsl-1.15/

.PHONY: all clean


