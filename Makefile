CC=g++
CFLAGS=-Wall -g
IFLAGS=
LFLAGS=

all: branched

branched.o : branched.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c branched.cc

surface.o : surface.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c surface.cc

branched : branched.o surface.o
	$(CC) -o branched branched.o surface.o $(LFLAGS)

clean : 
	rm *.o
	rm branched
