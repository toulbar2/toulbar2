SRCS = $(wildcard tb2*.cpp)

OBJS = $(SRCS:.cpp=.o)

LIBSRCS1 = $(SRCS:tb2main.cpp=)
LIBSRCS = $(LIBSRCS1:tb2solver.cpp=)

LIBOBJS = $(LIBSRCS:.cpp=.o)

CC = g++
ICC = icc

# Library compatible with Ilog Solver:
# In order to link ToulBar2 with Ilog Solver, uncomment the following two lines starting with ILOGFLAGS and ILOGLIBS
# Then recompile libtb2.so (make clean;make libtb2.so) and finaly compile iloglue (make iloglue)
# Usage: iloglue problem.wcsp [verbosity]
#ILOGFLAGS = -DILOGLUE -DIL_STD -I/usr/local/Ilog/solver60/include -I/usr/local/Ilog/concert20/include
#ILOGLIBS = -L/usr/local/Ilog/solver60/lib/i86_linux2_glibc2.3_gcc3.2/static_pic -L/usr/local/Ilog/concert20/lib/i86_linux2_glibc2.3_gcc3.2/static_pic -lsolverfloat -lsolver -lconcert

CFLAGS = -Wall $(ILOGFLAGS) -DWCSPFORMATONLY
# -g3 -ftrapv

# -ftrapv -g3 -ansi -pedantic
#OFLAGS = -O2 -DNDEBUG
OFLAGS = -O3 -funroll-loops -fomit-frame-pointer -DNDEBUG
ICCFLAGS = -O3 -ipo -static -DNDEBUG -pg

PFLAGS = 
#-pg -g3

BUILD = $(CC) -shared

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) -c -o $@  $< $(PFLAGS)

toulbar2: $(OBJS)
	$(CC) $(OBJS) -o $@ $(PFLAGS)

toulbar2i: $(SRCS)
	$(ICC) $(SRCS) -o $@ $(ICCFLAGS) $(CFLAGS)

libtb2.so: $(LIBOBJS)
	$(BUILD) $(LIBOBJS) -o $@

iloglue: iloglue.cpp libtb2.so
	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< $(ILOGLIBS) -L. -ltb2

clean:
	-rm -f *.o toulbar2 toulbar2i libtb2.so iloglue

all: toulbar2
