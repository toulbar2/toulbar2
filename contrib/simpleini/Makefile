
CC = g++

CFLAGS = -std=c++17 -g -Wall -I../../src -DBOOST -DLONGLONG_COST -DLONGDOUBLE_PROB -DWCSPFORMATONLY

OFLAGS = -O3 -DNDEBUG

toulbar2test: toulbar2test.cpp ../../lib/Linux/libtb2.so
	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< -L../.. -L../../lib/Linux -rdynamic ../../lib/Linux/libtb2.so -lgmp -lboost_graph -lboost_iostreams -lz -llzma -Wl,-rpath,../..:../../lib/Linux

clean:
	-rm -f *.o toulbar2test

all: toulbar2test
