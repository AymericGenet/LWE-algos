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


void bkw_mem_create(bkw_mem_t * bkw, lwe_t lwe, int a, int d, int m) {
    size_t i;
    char path[128];

    bkw->lwe = lwe;
    bkw->a = a;
    bkw->b = (lwe.n/((double) a) + 0.5); /* cheap round */
    bkw->d = d;
    bkw->m = m;

    bkw->tables = malloc(a * sizeof(FILE *));
    bkw->sample_pos = calloc(a, sizeof(vec_t));
    bkw->sample_neg = calloc(a, sizeof(vec_t));

    for (i = 0; i < a; ++i) {
        sprintf(path, "tables/lf2.T.%i.txt", i);
        open_table(&(bkw->tables[i]), path);
        bkw->sample_pos[i] = NULL;
        bkw->sample_neg[i] = NULL;
    }
}

void bkw_mem_free(bkw_mem_t * bkw) {
    size_t i;

    for (i = 0; i < bkw->a; ++i) {
        close_table(&(bkw->tables[i]));
        if (bkw->sample_pos[i] != NULL) {
            free(bkw->sample_pos[i]);
            free(bkw->sample_neg[i]);
        }
    }
    free(bkw->sample_pos);
    free(bkw->sample_neg);
    free(bkw->tables);
}

int bkw_mem_lf1(vec_t res, bkw_mem_t bkw, int l, math_t ** aux) {
    size_t i;
    int b, n; /*, d; */
    long q;
    FILE ** tables;

    /* for visibility concerns */
    n = bkw.lwe.n;
    q = bkw.lwe.q;
    b = bkw.b;
    /* d = bkw.d; */ /* TODO : d reduction */
    tables = bkw.tables;

    while (1) {
        /* samples from previous layer */
        if (l == 0) {
            return lwe_oracle(res, bkw.lwe);
        } else {
            if (!bkw_mem_lf1(aux[0], bkw, l - 1, aux + 3)) {
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
        rewind(tables[l - 1]);
        while (read_sample(aux[2], tables[l - 1], q, n + 1)) {
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
        append_sample(tables[l - 1], aux[0], n + 1);
    }
    return 0; /* should not happen */
}

int bkw_mem_lf2(vec_t res, bkw_mem_t bkw, int l, math_t ** aux) {
    size_t i;
    int n, b; /*, d; */
    long q;
    vec_t sample_pos;
    vec_t sample_neg;
    FILE ** tables;

    /* if call to lwe_oracle */
    if (l == 0) {
        return lwe_oracle(res, bkw.lwe);
    }

    /* for visibility concerns */
    n = bkw.lwe.n;
    q = bkw.lwe.q;
    b = bkw.b;
    /* d = bkw.d; */ /* TODO : d reduction */
    tables = bkw.tables;

    /* finds layer at which the last sample was output */
    sample_pos = bkw.sample_pos[l-1];
    sample_neg = bkw.sample_neg[l-1];

    /* upon first layer call */
    if (sample_pos == NULL || sample_neg == NULL) {
        bkw.sample_pos[l-1] = calloc(n + 1, sizeof(math_t));
        bkw.sample_neg[l-1] = calloc(n + 1, sizeof(math_t));
        sample_pos = bkw.sample_pos[l-1];
        sample_neg = bkw.sample_neg[l-1];
    } else {
        /* tries to find a further collision in the table */
        while (read_sample(aux[0], tables[l-1], q, n)) {
            if (equals(aux[0], sample_pos, (l-1) * b, l * b)) {
                for (i = 0; i < n + 1; ++i) {
                    res[i] = (sample_neg[i] + aux[0][i]) % q;
                }
                return 1; /* true */
            }
            if (equals(aux[0], sample_neg, (l-1) * b, l * b)) {
                for (i = 0; i < n + 1; ++i) {
                    res[i] = (aux[0][i] + sample_pos[i]) % q;
                }
                return 1; /* true */
            }
        }

        /* otherwise, stores it */
        append_sample(tables[l-1], sample_pos, n + 1);
    }

    while (1) {
        /* samples from previous layer */
        if (!bkw_mem_lf2(sample_pos, bkw, l - 1, aux + 1)) {
            return 0; /* false */
        }

        /* negates sample */
        for (i = 0; i < n + 1; ++i) {
            sample_neg[i] = (q - sample_pos[i]) % q;
        }

        /* resets reading pointer */
        rewind(tables[l-1]);

        /* if first elements already 0, returns it */
        if (zero(sample_pos, (l-1) * b, l * b)) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = sample_pos[i];
            }

            return 1; /* true */
        }

        /* checks for collision for either a or -a */
        while (read_sample(aux[0], tables[l-1], q, n)) {
            if (equals(aux[0], sample_pos, (l-1) * b, l * b)) {
                for (i = 0; i < n + 1; ++i) {
                    res[i] = (sample_neg[i] + aux[0][i]) % q;
                }
                return 1; /* true */
            }
            if (equals(aux[0], sample_neg, (l-1) * b, l * b)) {
                for (i = 0; i < n + 1; ++i) {
                    res[i] = (aux[0][i] + sample_pos[i]) % q;
                }
                return 1; /* true */
            }
        }

        /* otherwise, stores it, and repeat */
        append_sample(tables[l-1], sample_pos, n + 1);
    }

    return 0;
}
