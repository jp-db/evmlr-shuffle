CC = gcc
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
LIBS = -lflint -lgmp
UTILS = evmlr_utils.c
BENCH = bench.c cpucycles.c
TEST = test.c

evmlr_hpke_EXTRA = evmlr_otse.o evmlr_mlpke.o

.PHONY: all clean encrypt commit

all: mlpke commit otse hpke

mlpke: evmlr_mlpke.bin
	./$<

commit: evmlr_commit.bin
	./$<

otse: evmlr_otse.bin
	./$<

hpke: $(evmlr_hpke_EXTRA) evmlr_hpke.bin
	./$(word $(words $^), $^)

%.bin: %.c $(UTILS) $(TEST) $(BENCH)
	$(CC) $(CFLAGS) -DMAIN $^ $($*_EXTRA) -o $@ $(LIBS)
	chmod +x $@ # make executable

%.o: %.c $(UTILS)
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

clean:
	rm -f *.o *.so *.bin