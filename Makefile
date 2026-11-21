# Build knobs. Override on the command line, e.g.
#   make ARCH=            # portable binaries (no -march=native)
#   make OPT=-O0          # unoptimised, for debugging
#   make EXTRA_CFLAGS=-fsanitize=address
ARCH        ?= -march=native -mtune=native
OPT         ?= -O3
EXTRA_CFLAGS ?=

CC  = gcc
CXX = g++

CFLAGS   = $(OPT) $(ARCH) -Wall -ggdb -pthread $(EXTRA_CFLAGS)
DEPFLAGS = -MMD -MP
LIBS     = -lflint -lgmp -lstdc++ -lm

PREFIX  = evmlr
TARGETS = mlpke commit otse hpke shuffle lin_proof lin_comp bin_proof voting enc_proof serialize
BINS    = $(TARGETS:%=$(PREFIX)_%.bin)

# Objects every binary needs: shared helpers plus the test/benchmark harness.
COMMON_OBJS = $(PREFIX)_utils.o $(PREFIX)_crt.o test.o bench.o cpucycles.o randombytes.o

# Extra objects required by individual binaries.
$(PREFIX)_commit_EXTRA    = fastrandombytes.o
$(PREFIX)_hpke_EXTRA      = $(PREFIX)_otse.o $(PREFIX)_mlpke.o
$(PREFIX)_shuffle_EXTRA   = $(PREFIX)_hpke.o $(PREFIX)_commit.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o sha224-256.o fastrandombytes.o $(PREFIX)_lin_proof.o $(PREFIX)_bin_proof.o $(PREFIX)_challenge.o gaussian.o
$(PREFIX)_lin_proof_EXTRA = sha224-256.o fastrandombytes.o $(PREFIX)_challenge.o gaussian.o
$(PREFIX)_lin_comp_EXTRA  = sha224-256.o fastrandombytes.o $(PREFIX)_challenge.o
$(PREFIX)_bin_proof_EXTRA = sha224-256.o fastrandombytes.o $(PREFIX)_commit.o $(PREFIX)_challenge.o gaussian.o
$(PREFIX)_voting_EXTRA    = $(PREFIX)_shuffle.o $(PREFIX)_hpke.o $(PREFIX)_commit.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o sha224-256.o fastrandombytes.o $(PREFIX)_lin_proof.o $(PREFIX)_bin_proof.o $(PREFIX)_challenge.o $(PREFIX)_enc_proof.o gaussian.o
$(PREFIX)_enc_proof_EXTRA = $(PREFIX)_lin_proof.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o $(PREFIX)_hpke.o $(PREFIX)_challenge.o sha224-256.o fastrandombytes.o gaussian.o
$(PREFIX)_serialize_EXTRA = $(PREFIX)_voting.o $(PREFIX)_shuffle.o $(PREFIX)_hpke.o $(PREFIX)_commit.o $(PREFIX)_mlpke.o $(PREFIX)_otse.o sha224-256.o fastrandombytes.o $(PREFIX)_lin_proof.o $(PREFIX)_bin_proof.o $(PREFIX)_challenge.o $(PREFIX)_enc_proof.o gaussian.o

.PHONY: all build test bench run clean help shared $(TARGETS)

# Shared library: the protocol without the test/benchmark harness. Its objects
# are compiled again with -fPIC (as %.pic.o) so the binaries stay as they are.
LIB_NAME    = lib$(PREFIX).so
LIB_VERSION = 0
LIB_SONAME  = $(LIB_NAME).$(LIB_VERSION)
LIB_OBJS    = $(addsuffix .pic.o, \
	$(PREFIX)_utils $(PREFIX)_crt $(PREFIX)_mlpke $(PREFIX)_otse $(PREFIX)_hpke \
	$(PREFIX)_commit $(PREFIX)_challenge $(PREFIX)_lin_proof $(PREFIX)_bin_proof \
	$(PREFIX)_enc_proof $(PREFIX)_shuffle $(PREFIX)_voting $(PREFIX)_serialize \
	sha224-256 fastrandombytes gaussian randombytes cpucycles)

## all:   build every binary (does not run anything)
all: build

## build: compile every binary
build: $(BINS)

## test:  build, then run the test suites only (fast)
test: $(BINS)
	@fail=0; for b in $(BINS); do \
		echo "=== $$b ==="; ./$$b test || fail=1; \
	done; exit $$fail

## bench: build, then run the benchmark suites only (slow)
bench: $(BINS)
	@fail=0; for b in $(BINS); do \
		echo "=== $$b ==="; ./$$b bench || fail=1; \
	done; exit $$fail

## shared: build the shared library lib$(PREFIX).so (link with -levmlr)
shared: $(LIB_NAME)

$(LIB_NAME): $(LIB_SONAME)
	ln -sf $< $@

$(LIB_SONAME): $(LIB_OBJS)
	$(CC) $(CFLAGS) -shared -Wl,-soname,$(LIB_SONAME) -Wl,--no-undefined $^ -o $@ $(LIBS)

## run:   build, then run tests and benchmarks for every binary
run: $(BINS)
	@fail=0; for b in $(BINS); do \
		echo "=== $$b ==="; ./$$b all || fail=1; \
	done; exit $$fail

## <name>: build and run one binary, e.g. `make shuffle` (see TARGETS below)
$(TARGETS): %: $(PREFIX)_%.bin
	./$<

# Some sources use x86 intrinsics unconditionally and therefore need those
# instruction sets enabled even when ARCH does not imply them (e.g. ARCH=).
fastrandombytes.o: EXTRA_CFLAGS += -maes -msse4.1
gaussian.o:        EXTRA_CFLAGS += -msse4.1
fastrandombytes.pic.o: EXTRA_CFLAGS += -maes -msse4.1
gaussian.pic.o:        EXTRA_CFLAGS += -msse4.1

.SECONDARY: # keep the %.main.o intermediates so rebuilds stay incremental

.SECONDEXPANSION: # expand $(<name>_EXTRA) after the stem is known
%.bin: %.main.o $$($$*_EXTRA) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

# The translation unit that provides main(); compiled with -DMAIN so that the
# test/benchmark section of the file is included.
%.main.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -DMAIN -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

%.pic.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -fPIC -c $< -o $@

%.pic.o: %.cpp
	$(CXX) $(CFLAGS) $(DEPFLAGS) -fPIC -c $< -o $@

help:
	@echo "Targets:"
	@sed -n 's/^## //p' $(MAKEFILE_LIST)
	@echo
	@echo "TARGETS = $(TARGETS)"

clean:
	rm -f *.o *.d *.so *.so.* *.bin

# Header dependencies recorded by -MMD, so that editing e.g. evmlr_params.h
# rebuilds everything that includes it.
-include $(wildcard *.d)
