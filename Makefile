CC=g++
CFLAGS=-O3 #-Wall -g
GFLAGS=-I/usr/X11R6/include
IFLAGS=-I/sw/include
LFLAGS=-L/usr/X11R6/lib -lX11 -L/sw/lib -lgmp

all: branched

branched.o : branched.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c branched.cc

surface.o : surface.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c surface.cc

graphics.o : graphics.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graphics.cc

rational.o : rational.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c rational.cc

perm.o : perm.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c perm.cc

branched : branched.o surface.o graphics.o rational.o perm.o
	$(CC) -o branched branched.o surface.o graphics.o rational.o perm.o $(LFLAGS)

clean : 
	rm *.o
	rm branched
