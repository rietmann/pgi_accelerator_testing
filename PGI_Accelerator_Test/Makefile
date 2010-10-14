#
# Makefile for Patus stencil benchmark
#
# Note: $(PATUS_*) variables will be automatically replaced by the
# required runtime files by Patus.
#

CC = pgcc
CFLAGS = -O3 -ta=nvidia -Minfo
CFLAGS_CPU = -O3 -Minfo

bench: kernel.o driver.o timer.o
	$(CC) $(CFLAGS) -o $@ $+

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *.o bench
