# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
CFLAGS = -o
LIBS = -lm

# Output directory
PROGS = multipolator

# Name of the executables:
all: $(PROGS)

multipolator: ATMO_interpolate.c multipolator_w_printing.c multipolator.h
	$(CC) $(CFLAGS) -lm -o $@ ATMO_interpolate.c multipolator_w_printing.c

clean:
	rm -f $(PROGS)
