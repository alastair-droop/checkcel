CC=gcc
CFLAGS=-Wall

all:
	$(CC) $(CFLAGS) -o checkcel *.c
	
clean:
	-rm checkcel
