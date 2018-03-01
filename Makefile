CC=icc
CFLAGS=-std=gnu99 -fopenmp
LDFLAGS=
LDLIBS=-lm -lfftw3 -lqdsp

.PHONY: all
all: oldpic ftpic

.PHONY: debug
debug: oldpic ftpic
debug: CFLAGS += -g -O0

oldpic: oldpic.c
	$(CC) -o oldpic $(CFLAGS) $(LDFLAGS) oldpic.c $(LDLIBS)

ftpic: ftpic.c
	$(CC) -o ftpic $(CFLAGS) $(LDFLAGS) ftpic.c libusfft.a -lgfortran $(LDLIBS)
