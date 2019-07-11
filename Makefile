# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
CFLAGS = -o
LIBS = -lm

# Output directory
PROGS = multipolator

# Name of the executables:
all: $(PROGS)

multipolator: ATMO_interpolate.c multipolator.c multipolator.h
	$(CC) $(CFLAGS) -lm -o $@ ATMO_interpolate.c multipolator.c

clean:
	rm -f $(PROGS)
