SRCSFULL = $(wildcard tb2*.cpp applis/*.cpp core/*.cpp globals/*.cpp incop/*.cpp search/*.cpp utils/*.cpp vns/*.cpp ToulbarVersion.cpp)

# Comment the following line if BOOST Graph library uninstalled
BOOST = -DBOOST

# Same thing for XML
#XML = XMLFLAG
#XMLI = -I/usr/include/libxml2 -DXMLFLAG
# -DXMLFLAG -DMAXCSP

ifeq ($(BOOST),-DBOOST)
SRCS = $(SRCSFULL)
else
SRCS = $(SRCSFULL:utils/tb2boostgraph.cpp=)
endif

ifeq ($(XML),XMLFLAG)
OBJS = $(SRCS:.cpp=.o) -lxml2 -lboost_graph -lboost_iostreams -lgmp -lz -llzma
# -lxml2 -lz -llzma -lm -lpthread
else
OBJS = $(SRCS:.cpp=.o) -lboost_graph -lboost_iostreams -lgmp -lz -llzma
endif

LIBSRCS = $(SRCS:tb2main.cpp=)
LIBOBJS = $(LIBSRCS:.cpp=.o)

CC = g++

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

## PLPS cost functions needs CPLEX installed (update the following path depending on your configuration and uncomment the following three lines) *** Warning! Beta untested version for PLPS ***
##CPLEX_LOCATION = /opt/ibm/ILOG/CPLEX_Studio1210
##CPLEXFLAGS = -DILOGCPLEX -DIL_STD -fPIC -fexceptions -fno-strict-aliasing -pthread -I$(CPLEX_LOCATION)/cplex/include/ilcplex -I$(CPLEX_LOCATION)/cplex/include -I$(CPLEX_LOCATION)/concert/include
##CPLEXLIB =  -L$(CPLEX_LOCATION)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX_LOCATION)/concert/lib/x86-64_linux/static_pic -lilocplex -lcplex -lconcert -pthread -rdynamic 

CFLAGS = -std=c++11 -g -Wall -I. $(ILOGFLAGS) $(BOOST) $(XMLI) -DLONGLONG_COST -DLONGDOUBLE_PROB -DWCSPFORMATONLY
# -DLONGLONG_COST # uses type long long for encoding costs, to be replaced by -DINT_COST if memory is an issue (but costs are limited to 119304647)
# -DLONGDOUBLE_PROB # uses type long double for encoding (log) probabilities instead of type double if using -DDOUBLE_PROB
# -DWCSPFORMATONLY # assumes positive domains starting at value 0 (as in the wcsp format file) *** Incompatible with Numberjack ***
# -DNUMBERJACK # 
# -DOPENMPI # enables variable neighborhood search
# -fsanitize=signed-integer-overflow -fsanitize=float-cast-overflow -g3 # catches numerical integer overflows
# -Wconversion # catch type conversion from larger to smaller type (e.g., Cost to int)

OFLAGS = -O3 -DNDEBUG
# -O3 -DNDEBUG

PFLAGS =
# -pg -g3

BUILD = $(CC) -shared

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp %.hpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp %.h
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

toulbar2: $(OBJS)
	$(CC) $(OBJS) -o $@ $(PFLAGS) -static -static-libgcc

libtb2.so: LIBFLAGS = -fPIC
libtb2.so: $(LIBOBJS)
	$(BUILD) $(LIBFLAGS) $(LIBOBJS) -o $@

toulbar2test: toulbar2test.cpp libtb2.so
	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< ./libtb2.so -lboost_graph -lboost_iostreams -lgmp -lz -llzma

iloglue: iloglue.cpp libtb2.so
	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< $(ILOGLIBS) -L. -ltb2 -lboost_graph -lboost_iostreams -lgmp -lz -llzma

clean:
	-rm -f *.o applis/*.o core/*.o globals/*.o incop/*.o search/*.o utils/*.o vns/*.o toulbar2 libtb2.so iloglue toulbar2test

all: toulbar2
