CC=gcc
CFLAGS=-O3

all: main

main: src/main.c
	$(CC) $(LIB) $(CFLAGS) -fopenmp $< -o $@
