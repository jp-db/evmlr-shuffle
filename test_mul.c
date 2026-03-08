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

    nmod_poly_t A, B, C_flint, C_crt;
    nmod_poly_init(A, 3109);
    nmod_poly_init(B, 3109);
    nmod_poly_init(C_flint, 3109);
    nmod_poly_init(C_crt, 3109);

    flint_rand_t state;
    flint_rand_init(state);
    
    evmlr_poly_crt_t crt_A, crt_B, crt_C;
    evmlr_poly_crt_init(crt_A, 3109);
    evmlr_poly_crt_init(crt_B, 3109);
    evmlr_poly_crt_init(crt_C, 3109);

    for (int i = 0; i < 10000; i++) {
        nmod_poly_randtest(A, state, 512);
        nmod_poly_randtest(B, state, 512);
        
        nmod_poly_mulmod(C_flint, A, B, cyclo);
        
        evmlr_poly_to_crt(crt_A, A, ctx);
        evmlr_poly_to_crt(crt_B, B, ctx);
        evmlr_poly_crt_mul(crt_C, crt_A, crt_B, ctx);
        evmlr_poly_from_crt(C_crt, crt_C, ctx);
        
        if (!nmod_poly_equal(C_flint, C_crt)) {
            printf("Mismatch on iter %d!\n", i);
            return 1;
        }
    }
    
    printf("10000 CRTs multiplications match perfectly.\n");
    return 0;
}
