CC = gcc
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
LIBS = -lflint -lgmp
UTILS = evmlr_utils.c
BENCH = bench.c cpucycles.c
TEST = test.c

.PHONY: all clean encrypt commit

all: encrypt commit

encrypt: evmlr_enc.o
	./evmlr_enc.o

commit: evmlr_commit.o
	./evmlr_commit.o

%.o: %.c $(UTILS) $(TEST) $(BENCH)
	$(CC) $(CFLAGS) -DMAIN $^ -o $@ $(LIBS)
	chmod +x $@ # make executable

clean:
	rm -f *.o