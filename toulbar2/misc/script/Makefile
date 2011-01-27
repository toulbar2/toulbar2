SRCSFULL = $(wildcard tb2*.cpp)

# Remove the first character in the following line if compile with BOOST Graph library
BOOST = -DLINUX -DLONGLONG_COST -DWIDE_STRING -DLONGDOUBLE_PROB

# Same thing for XML
#XML = XMLFLAG
#XMLI = -I/usr/include/libxml2 -DXMLFLAG
# -DXMLFLAG -DMAXCSP

ifeq ($(BOOST),-DBOOST)
SRCS = $(SRCSFULL)
else
SRCS = $(SRCSFULL:tb2boostgraph.cpp=)
endif

ifeq ($(XML),XMLFLAG)
OBJS = $(SRCS:.cpp=.o) -lxml2
# -lxml2 -lz -lm -lpthread
else
OBJS = $(SRCS:.cpp=.o) -lgmp
endif

LIBSRCS1 = $(SRCS:tb2main.cpp=)
LIBSRCS2 = $(LIBSRCS1:tb2solver.cpp=)
LIBSRCS = $(LIBSRCS2:tb2btd.cpp=)
LIBOBJS = $(LIBSRCS:.cpp=.o)

CC = g++
ICC = icc

# Library compatible with Ilog Solver:
# In order to link ToulBar2 with Ilog Solver, uncomment the following two lines starting with ILOGFLAGS and ILOGLIBS, and uncomment #define INT_COST (also adding comment to #define LONGLONG_COST) in tb2types.hpp
# Then recompile libtb2.so (make clean;make libtb2.so) and finaly compile iloglue (make iloglue)
# Usage: iloglue 16wqueens.wcsp [verbosity level]
# Ilog Solver 6.0
#ILOGFLAGS = -DILOGLUE -DIL_STD -I/usr/local/Ilog/solver60/include -I/usr/local/Ilog/concert20/include
#ILOGLIBS = -L/usr/local/Ilog/solver60/lib/i86_linux2_glibc2.3_gcc3.2/static_pic -L/usr/local/Ilog/concert20/lib/i86_linux2_glibc2.3_gcc3.2/static_pic -lsolverfloat -lsolver -lconcert
# Ilog Solver 6.5
# (warning: compile libtb2.so and iloglue with -O optimization and not -O3)
#ILOGFLAGS = -DILOGLUE -DIL_STD -I/usr/local/ILOG/solver65/include -I/usr/local/ILOG/concert25/include
#ILOGLIBS = -L/usr/local/ILOG/solver65/lib/x86-64_rhel4.0_3.4/static_pic -L/usr/local/ILOG/concert25/lib/x86-64_rhel4.0_3.4/static_pic -lsolverfloat -lsolver -lconcert

CFLAGS = -g -Wall $(ILOGFLAGS) -DWCSPFORMATONLY -DNARYCHAR $(BOOST) $(XMLI)
# -O1 -g3 -ftrapv -ansi -pedantic
# -O1 -g3 -ftrapv -ansi -pedantic -Wall -Wno-long-long
# -DNARYCHAR  # use this option to print nary constraints (including separators) correctly

OFLAGS = -O3 -DNDEBUG
# -O1 -DNDEBUG
# -O2 -DNDEBUG

ICCFLAGS = 
#-O3 -ipo -static -DNDEBUG

PFLAGS = -pg -g3

BUILD = $(CC) -shared

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp %.hpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

toulbar2: $(OBJS)
	$(CC) $(OBJS) -o $@ $(PFLAGS)
#	$(CC) $(OBJS) -o $@ $(PFLAGS) -static -static-libgcc

toulbar2i: $(SRCS)
	$(ICC) $(SRCS) -o $@ $(ICCFLAGS) $(CFLAGS)

libtb2.so: LIBFLAGS = -fPIC
libtb2.so: $(LIBOBJS)
	$(BUILD) $(LIBFLAGS) $(LIBOBJS) -o $@

iloglue: iloglue.cpp libtb2.so
	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< $(ILOGLIBS) -L. -ltb2

clean:
	-rm -f *.o toulbar2 toulbar2i libtb2.so iloglue

all: toulbar2
