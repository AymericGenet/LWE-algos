/*
 * bkw_mem.c
 *
 *  Created on: Apr 28, 2016
 *      Author: Aymeric Genet
 */

#include "bkw_mem.h"
#include "math.h"
#include "lwe.h"
#include "misc.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int bkw_mem_lf1(vec_t res, int n, long q, int b, int l, FILE ** table,
                math_t ** aux) {
    size_t i;

    while (1) {
        /* samples from previous layer */
        if (l == 0) {
            return lwe_oracle(res, n, q);
        } else {
            if (!bkw_mem_lf1(aux[0], n, q, b, l - 1, table, aux + 3)) {
                return 0; /* false */
            }
        }

        /* if first elements already 0, returns it */
        if (zero(aux[0], (l - 1) * b, l * b)) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = aux[0][i];
            }
            return 1; /* true */
        }

        /* negates sample */
        for (i = 0; i < n + 1; ++i) {
            aux[1][i] = (q - aux[0][i]) % q;
        }

        /* checks for collision for either a or -a */
        rewind(table[l - 1]);
        while (read_sample(aux[2], table[l - 1], q, n + 1)) {
            if (equals(aux[0], aux[2], (l-1) * b, l * b)) {
                for (i = 0; i < n + 1; ++i) {
                    res[i] = (aux[1][i] + aux[2][i]) % q;
                }
                return 1; /* true */
            }
            if (equals(aux[1], aux[2], (l-1) * b, l * b)) {
                for (i = 0; i < n + 1; ++i) {
                    res[i] = (aux[0][i] + aux[2][i]) % q;
                }
                return 1; /* true */
            }
        }

        /* otherwise, stores it, and repeat */
        append_sample(table[l - 1], aux[0], n + 1);
    }
    return 0; /* should not happen */
}

int bkw_mem_lf2(vec_t res, int l, bkw_mem_t * bkw, math_t ** aux) {
    size_t i;
    vec_t sample_pos;
    vec_t sample_neg;

    /* if call to lwe_oracle */
    if (l == 0) {
        return lwe_oracle(res, bkw->n, bkw->q);
    }

    /* finds layer at which the last sample was output */
    sample_pos = bkw->sample_pos[l-1];
    sample_neg = bkw->sample_neg[l-1];

    /* upon first layer call */
    if (sample_pos == NULL || sample_neg == NULL) {
        bkw->sample_pos[l-1] = calloc(bkw->n + 1, sizeof(math_t));
        bkw->sample_neg[l-1] = calloc(bkw->n + 1, sizeof(math_t));
        sample_pos = bkw->sample_pos[l-1];
        sample_neg = bkw->sample_neg[l-1];
    } else {
        /* tries to find a further collision in the table */
        while (read_sample(aux[0], bkw->tables[l-1], bkw->q, bkw->n)) {
            if (equals(aux[0], sample_pos, (l-1) * bkw->b, l * bkw->b)) {
                for (i = 0; i < bkw->n + 1; ++i) {
                    res[i] = (sample_neg[i] + aux[0][i]) % bkw->q;
                }
                return 1; /* true */
            }
            if (equals(aux[0], sample_neg, (l-1) * bkw->b, l * bkw->b)) {
                for (i = 0; i < bkw->n + 1; ++i) {
                    res[i] = (aux[0][i] + sample_pos[i]) % bkw->q;
                }
                return 1; /* true */
            }
        }

        /* otherwise, stores it */
        append_sample(bkw->tables[l-1], sample_pos, bkw->n + 1);
    }

    while (1) {
        /* samples from previous layer */
        if (!bkw_mem_lf2(sample_pos, l - 1, bkw, aux + 1)) {
            return 0; /* false */
        }

        /* negates sample */
        for (i = 0; i < bkw->n + 1; ++i) {
            sample_neg[i] = (bkw->q - sample_pos[i]) % bkw->q;
        }

        /* resets reading pointer */
        rewind(bkw->tables[l-1]);

        /* if first elements already 0, returns it */
        if (zero(sample_pos, (l-1) * bkw->b, l * bkw->b)) {
            for (i = 0; i < bkw->n + 1; ++i) {
                res[i] = sample_pos[i];
            }

            return 1; /* true */
        }

        /* checks for collision for either a or -a */
        while (read_sample(aux[0], bkw->tables[l-1], bkw->q, bkw->n)) {
            if (equals(aux[0], sample_pos, (l-1) * bkw->b, l * bkw->b)) {
                for (i = 0; i < bkw->n + 1; ++i) {
                    res[i] = (sample_neg[i] + aux[0][i]) % bkw->q;
                }
                return 1; /* true */
            }
            if (equals(aux[0], sample_neg, (l-1) * bkw->b, l * bkw->b)) {
                for (i = 0; i < bkw->n + 1; ++i) {
                    res[i] = (aux[0][i] + sample_pos[i]) % bkw->q;
                }
                return 1; /* true */
            }
        }

        /* otherwise, stores it, and repeat */
        append_sample(bkw->tables[l-1], sample_pos, bkw->n + 1);
    }

    return 0;
}
