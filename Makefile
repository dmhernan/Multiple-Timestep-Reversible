# Makefile for MTR (Multiple Timestep Reversible) integrator

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm

SRCS    = main.c integrator.c maps.c kepler.c coordinates.c physics.c init.c
OBJS    = $(SRCS:.c=.o)
TARGET  = mtr

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c mtr.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)
