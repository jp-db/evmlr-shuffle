#include "evmlr_params.h"
#include "flint/nf.h"
#include "flint/nf_elem.h"

// If we have lists of elements (a_1, ..., a_N ) and (b_1, ..., b_N ), where all the a_i, b_i
// are in some field F, then the two lists are a permutation of each other if and only
// if over F[X], we have:
// \sum_{i=1}^{N} (a_i - X) = \sum_{i=1}^{N} (b_i - X)

int is_perm_over_f(nf_elem_t *a, nf_elem_t *b, slong len, nf_t field) {
    nf_elem_t sum_a, sum_b;
    nf_elem_init(sum_a, field);
    nf_elem_init(sum_b, field);
    nf_elem_zero(sum_a);
    nf_elem_zero(sum_b);

    nf_elem_t temp_a, temp_b;
    nf_elem_init(temp_a, field);
    nf_elem_init(temp_b, field);
    for (slong i = 0; i < len; i++) {
        // TODO
    }

    return nf_elem_equal(sum_a, sum_b, field);
}