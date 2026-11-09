/* Adapted from NewHope 
 * https://github.com/tpoeppelmann/newhope */

#include "randombytes.h"

#ifdef __EMSCRIPTEN__
#include <stdio.h>
#include <stdlib.h>
#include <sys/random.h>

/* WebAssembly has no /dev/urandom unless the virtual filesystem is linked in;
 * getentropy() reaches crypto.getRandomValues() (browser) or
 * crypto.randomFillSync() (Node) directly, at most 256 bytes per call. */
int randombytes(unsigned char *x, unsigned long long xlen) {
	while (xlen > 0) {
		size_t n = xlen < 256 ? (size_t) xlen : 256;

		if (getentropy(x, n) != 0) {
			perror("randombytes: getentropy");
			abort();
		}
		x += n;
		xlen -= n;
	}

	return 0;
}

#else
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

static int fd = -1;

int randombytes(unsigned char *x, unsigned long long xlen) {
	int i;

	if (fd == -1) {
		for (;;) {
			fd = open("/dev/urandom", O_RDONLY);
			if (fd != -1)
				break;
			sleep(1);
		}
	}

	while (xlen > 0) {
		if (xlen < 1048576)
			i = xlen;
		else
			i = 1048576;

		i = read(fd, x, i);
		if (i < 1) {
			sleep(1);
			continue;
		}

		x += i;
		xlen -= i;
	}

	return 0;
}
#endif
