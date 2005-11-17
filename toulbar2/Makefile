SRCS = $(wildcard toulbar2/tb2*.cpp)

OBJS = $(SRCS:.cpp=.o)

CC = g++
ICC = icc

CFLAGS = -Wall

# -ftrapv -g3 -ansi
OFLAGS = -O3 -funroll-loops -fomit-frame-pointer -DNDEBUG
ICCFLAGS = -O3 -ipo -static -DNDEBUG

PFLAGS = 
# -pg -g3

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) -c -o $@  $< $(PFLAGS)

toulbar2: $(OBJS)
	$(CC) $(OBJS) -o toulbar2/$@ $(PFLAGS)

toulbar2i: $(SRCS)
	$(ICC) $(SRCS) -o toulbar2/$@ $(ICCFLAGS) $(CFLAGS)

lint: $(SRCS)
	splint -nullassign -nullstate -globstate -compdef -unqualifiedtrans $(SRCS)

clean:
	-rm -f toulbar2/*.o toulbar2/toulbar2 toulbar2/toulbar2i

all:	toulbar2 toulbar2i
