CC=g++
CFLAGS=-Wall -g
GFLAGS=-I/usr/X11R6/include
IFLAGS=
LFLAGS=-L/usr/X11R6/lib -lX11

all: branched

branched.o : branched.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c branched.cc

surface.o : surface.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c surface.cc

graphics.o : graphics.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graphics.cc

branched : branched.o surface.o graphics.o
	$(CC) -o branched branched.o surface.o graphics.o $(LFLAGS)

clean : 
	rm *.o
	rm branched
