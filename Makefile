CC=gcc
CFLAGS=-std=gnu99 -fopenmp
LDFLAGS=
LDLIBS=-lm -lfftw3 -lqdsp

.PHONY: all
all: oldpic ftpic

.PHONY: debug
debug: oldpic ftpic
debug: CFLAGS += -g -O0

oldpic: oldpic.c
	gcc $(CFLAGS) $(LDFLAGS) -o oldpic oldpic.c $(LDLIBS)

ftpic: ftpic.c
	gcc $(CFLAGS) $(LDFLAGS) -o ftpic ftpic.c $(LDLIBS)
