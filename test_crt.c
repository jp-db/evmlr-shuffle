#include "evmlr_crt.h"
#include <stdio.h>
#include <assert.h>

int main() {
    evmlr_crt_ctx_t ctx;
    nmod_poly_t cyclo;
    nmod_poly_init(cyclo, 3109);
    nmod_poly_set_coeff_ui(cyclo, 512, 1);
    nmod_poly_set_coeff_ui(cyclo, 0, 1);
    
    evmlr_crt_setup(ctx, cyclo);

    nmod_poly_t A, B, C;
    nmod_poly_init(A, 3109);
    nmod_poly_init(B, 3109);
    nmod_poly_init(C, 3109);

    flint_rand_t state;
    flint_rand_init(state);
    
    for (int i = 0; i < 1000; i++) {
        nmod_poly_randtest(A, state, 512);
        
        evmlr_poly_crt_t crt_A;
        evmlr_poly_crt_init(crt_A, 3109);
        
        evmlr_poly_to_crt(crt_A, A, ctx);
        evmlr_poly_from_crt(B, crt_A, ctx);
        
        if (!nmod_poly_equal(A, B)) {
            printf("Mismatch!\n");
            return 1;
        }
        
        evmlr_poly_crt_clear(crt_A);
    }
    
    printf("1000 CRTs mapped back successfully.\n");
    return 0;
}
