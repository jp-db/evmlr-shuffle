CC = gcc
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
LIBS = -lflint -lgmp
PREFIX = evmlr
UTILS = $(PREFIX)_utils.c
BENCH = bench.c cpucycles.c
TEST = test.c
TARGETS = mlpke commit otse hpke shuffle lin_proof

$(PREFIX)_hpke_EXTRA = $(PREFIX)_otse.o $(PREFIX)_mlpke.o
$(PREFIX)_shuffle_EXTRA = $(PREFIX)_hpke.o $(PREFIX)_commit.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o sha224-256.o fastrandombytes.o $(PREFIX)_lin_proof.o
$(PREFIX)_lin_proof_EXTRA = sha224-256.o fastrandombytes.o

.PHONY: all clean $(TARGETS)

all: $(TARGETS)

$(TARGETS): %: $(PREFIX)_%.bin
	./$(word $(words $^), $^) # run last target (the binary)

.SECONDEXPANSION: # makefile magic to expand and get the extra dependencies
%.bin: %.c $$($$*_EXTRA) $(UTILS) $(TEST) $(BENCH)
	$(CC) $(CFLAGS) -DMAIN $^ -o $@ $(LIBS)
	chmod +x $@ # make executable

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

clean:
	rm -f *.o *.so *.bin