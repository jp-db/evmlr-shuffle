#include "cpucycles.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>

/* WebAssembly has no cycle counter; report nanoseconds of wall time instead
 * (emscripten_get_now() is performance.now(), in milliseconds). */
unsigned long long cpucycles() {
	return (unsigned long long) (emscripten_get_now() * 1e6);
}
#else
unsigned long long cpucycles() {
	unsigned long long result;
	asm volatile (".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax":"=a"
			(result)::"%rdx");
	return result;
}
#endif
