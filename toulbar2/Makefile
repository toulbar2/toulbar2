SRCS = $(wildcard tb2*.cpp)

OBJS = $(SRCS:.cpp=.o)

CC = g++
ICC = icc

CFLAGS = -Wall -g3

# -ftrapv -g3 -ansi -pedantic
OFLAGS = -O3 -funroll-loops -DNDEBUG
#OFLAGS = -O3 -funroll-loops -fomit-frame-pointer -DNDEBUG
ICCFLAGS = -O3 -ipo -static -DNDEBUG -pg

PFLAGS = -pg -g3

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) -c -o $@  $< $(PFLAGS)

toulbar2: $(OBJS)
	$(CC) $(OBJS) -o $@ $(PFLAGS)

toulbar2i: $(SRCS)
	$(ICC) $(SRCS) -o $@ $(ICCFLAGS) $(CFLAGS)

lint: $(SRCS)
	splint -nullassign -nullstate -globstate -compdef -unqualifiedtrans $(SRCS)

clean:
	-rm -f *.o toulbar2 toulbar2i

all:	toulbar2 toulbar2i
