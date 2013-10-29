CC=g++
CFLAGS=-wall -g
IFLAGS=
LFLAGS=

all: branched

branched.o : branched.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c branched.cc

branched : branched.o
	$(CC) -o branched branched.o $(LFLAGS)

clean : 
	rm *.o
	rm branched
