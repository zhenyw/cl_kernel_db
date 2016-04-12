
CFLAGS := -I/usr/include
LDFLAGS := -lOpenCL

all: run

run: run.c
	gcc -g -o run run.c $(CFLAGS) $(LDFLAGS)
