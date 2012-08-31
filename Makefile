CC = g++
CFLAGS = -Wall
AR = ar

all: util/single_1.so

util/single_1.so: src/single_1.pyx lib/quant.a setup.py
	python setup.py build_ext -i
	mv single_1.so util/

lib/quant.a: build/single_1.o
	$(AR) rcs lib/quant.a build/single_1.o

build/single_1.o: src/single_1.cc src/single_1.h
	$(CC) $(CFLAGS) -c src/single_1.cc -o build/single_1.o
	
clean:
	rm -rfv build/*
	rm -rfv util/single_1.so
	rm -rfv lib/quant.a
	rm -rfv src/single_1.cpp

.PHONY: all clean


