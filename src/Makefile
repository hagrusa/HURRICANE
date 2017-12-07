CC = gcc
CFLAGS = -g -std=c99 -O3
LDFLAGS = -lm

all: nbodymain

nbodymain: nbodymain.o bhtree.o integrators.o linkedlist.o
	$(CC) $(CFLAGS) -o nbodymain linkedlist.o bhtree.o integrators.o nbodymain.o $(LDFLAGS)

nbodymain.o: nbodymain.c
	$(CC) $(CFLAGS) -c nbodymain.c

bhtree.o: bhtree.c bhtree.h constants.h
	$(CC) $(CFLAGS) -c bhtree.c

linkedlist.o: linkedlist.c linkedlist.h
	$(CC) $(CFLAGS) -c linkedlist.c

integrators.o: integrators.c integrators.h constants.h
	$(CC) $(CFLAGS) -c integrators.c

clean:
	rm -f *.o nbodymain
