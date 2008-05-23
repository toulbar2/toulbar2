SRCSFULL = $(wildcard tb2*.cpp)

# Remove the first character in the following line if compile with BOOST Graph library
# BOOST = -DBOOST

# Same thing for XML
#XML = XMLFLAG
#XMLI = -I/usr/include/libxml2 -DXMLFLAG

ifeq ($(BOOST),-DBOOST)
SRCS = $(SRCSFULL)
else
SRCS = $(SRCSFULL:tb2boostgraph.cpp=)
endif

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

CFLAGS = -g -Wall $(ILOGFLAGS) -DWCSPFORMATONLY $(BOOST) $(XMLI)
# -g3 -ftrapv -ansi -pedantic
# -DNARYCHAR

OFLAGS = -O3 -DNDEBUG
# -O3 -DNDEBUG
# -O2

ICCFLAGS = 
#-O3 -ipo -static -DNDEBUG

PFLAGS = 
# -pg -g3
# -static -static-libgcc

ifeq ($(XML),XMLFLAG)
LIBFLAGS = -lxml2
else
LIBFLAGS = 
# -fPIC
endif

BUILD = $(CC) -shared

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp %.hpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

toulbar2: $(OBJS)
	$(CC) $(OBJS) -o $@ $(PFLAGS)

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
