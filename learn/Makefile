CC = g++

CFLAGS = -Wall -fPIC -Wconversion -Wextra -ggdb -DDEBUG \
-DHAVE_STD -DHAVE_NAMESPACES \
-Wall -W -Wmissing-prototypes \
-Wstrict-prototypes \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common

INCLUDES = -I../optpp-2.4/include -I../optpp-2.4/newmat11

LIBS = ../lib/libopt.a ../lib/libnewmat.a

all: optpp

optpp: optpp.cc
	$(CC) $(CFLAGS) $(INCLUDES) optpp.cc $(LIBS) -o optpp 
