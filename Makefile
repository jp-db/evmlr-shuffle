CC = gcc
CXX = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
LIBS = -lflint -lgmp -lstdc++ -lm
PREFIX = evmlr
UTILS = $(PREFIX)_utils.c
BENCH = bench.c cpucycles.c
TEST = test.c
TARGETS = mlpke commit otse hpke shuffle lin_proof lin_comp bin_proof voting enc_proof

$(PREFIX)_hpke_EXTRA = $(PREFIX)_otse.o $(PREFIX)_mlpke.o
$(PREFIX)_shuffle_EXTRA = $(PREFIX)_hpke.o $(PREFIX)_commit.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o sha224-256.o fastrandombytes.o $(PREFIX)_lin_proof.o $(PREFIX)_bin_proof.o $(PREFIX)_challenge.o gaussian.o
$(PREFIX)_lin_proof_EXTRA = sha224-256.o fastrandombytes.o $(PREFIX)_challenge.o gaussian.o
$(PREFIX)_lin_comp_EXTRA = sha224-256.o fastrandombytes.o $(PREFIX)_challenge.o
$(PREFIX)_bin_proof_EXTRA = sha224-256.o fastrandombytes.o $(PREFIX)_commit.o $(PREFIX)_challenge.o gaussian.o
$(PREFIX)_voting_EXTRA = $(PREFIX)_shuffle.o $(PREFIX)_hpke.o $(PREFIX)_commit.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o sha224-256.o fastrandombytes.o $(PREFIX)_lin_proof.o $(PREFIX)_bin_proof.o $(PREFIX)_challenge.o $(PREFIX)_enc_proof.o gaussian.o
$(PREFIX)_enc_proof_EXTRA = $(PREFIX)_lin_proof.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o $(PREFIX)_hpke.o $(PREFIX)_challenge.o sha224-256.o fastrandombytes.o gaussian.o

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

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@ $(LIBS)

clean:
	rm -f *.o *.so *.bin