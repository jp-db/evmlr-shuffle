CC = gcc
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
LIBS = -lflint -lgmp
UTILS = evmlr_utils.c

encrypt: evmlr_enc.c $(UTILS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
	chmod +x $@

clean:
	rm -f encrypt