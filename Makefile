CC=gcc
CFLAGS=-std=gnu99 -fopenmp
LDFLAGS=
LDLIBS=-lm -lfftw3 -lqdsp

.PHONY: all
all: oldpic ftpic oldpic2d ftpic2d

.PHONY: debug
debug: oldpic ftpic oldpic2d ftpic2d
debug: CFLAGS=-std=gnu99 -g -O0

oldpic: oldpic.o common.o
	$(CC) -o oldpic $(CFLAGS) $(LDFLAGS) oldpic.o common.o $(LDLIBS)

ftpic: ftpic.o common.o
	$(CC) -o ftpic $(CFLAGS) $(LDFLAGS) ftpic.o common.o libusfft.a -lgfortran $(LDLIBS)

oldpic2d: oldpic2d.o common2d.o
	$(CC) -o oldpic2d $(CFLAGS) $(LDFLAGS) oldpic2d.o common2d.o $(LDLIBS)

ftpic2d: ftpic2d.o common2d.o
	$(CC) -o ftpic2d $(CFLAGS) $(LDFLAGS) ftpic2d.o common2d.o libusfft.a -lgfortran $(LDLIBS)

common.o: common.h
oldpic.o: common.h
ftpic.o: common.h

common2d.o: common2d.h
oldpic2d.o: common2d.h
ftpic2d.o: common2d.h

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $<
