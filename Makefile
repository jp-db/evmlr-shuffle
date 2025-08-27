CC = gcc
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
LIBS = -lflint -lgmp
PREFIX = evmlr
UTILS = $(PREFIX)_utils.c
BENCH = bench.c cpucycles.c
TEST = test.c
TARGETS = mlpke commit otse hpke

$(PREFIX)_hpke_EXTRA = $(PREFIX)_otse.o $(PREFIX)_mlpke.o

.PHONY: all clean $(TARGETS)

all: $(TARGETS)

$(TARGETS): %: $(PREFIX)_%.bin
	./$(word $(words $^), $^) # run last target (the binary)

.SECONDEXPANSION: # makefile magic to expand and get the extra dependencies

%.bin: $$($$*_EXTRA) %.c $(UTILS) $(TEST) $(BENCH)
	$(CC) $(CFLAGS) -DMAIN $^ -o $@ $(LIBS)
	chmod +x $@ # make executable

%.o: $$($$*_EXTRA) %.c $(UTILS)
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

clean:
	rm -f *.o *.so *.bin