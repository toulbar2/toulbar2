CC=g++ -static
CFLAGS=-I. ./utils/
DEPS = utils/SimpleInit.h utils/SimpleGlob.h  utils/SimpleOpt.h utils/ConvertUTF.c

DEGUG = -DNDEBUG
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 
#	$(CC) -c -o $@ $< $(CFLAGS) -D DEBUG

cfntools: cfntools.o main_cfntools.o
	$(CC) -o cfntools cfntools.cpp main_cfntools.cpp 


.PHONY : clean
clean: 
	rm *.o $(objects)                              

