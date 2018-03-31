CC=icc
CFLAGS=-std=gnu99 -fopenmp
LDFLAGS=
LDLIBS=-lm -lfftw3 -lqdsp

.PHONY: all
all: oldpic ftpic

.PHONY: debug
debug: oldpic ftpic
debug: CFLAGS += -g -O0

oldpic: oldpic.o common.o
	$(CC) -o oldpic $(CFLAGS) $(LDFLAGS) oldpic.o common.o $(LDLIBS)

ftpic: ftpic.o common.o
	$(CC) -o ftpic $(CFLAGS) $(LDFLAGS) ftpic.o common.o libusfft.a -lgfortran $(LDLIBS)

common.o: common.h
oldpic.o: common.h
ftpic.o: common.h

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $<
