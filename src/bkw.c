/*
 * bkw.c
 *
 *  Created on: Apr 28, 2016
 *      Author: Aymeric Genet
 */

#include "bkw.h"
#include "math.h"
#include "lwe_oracle.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


static double * log_j = NULL;

int bkw_algo(math_t * res, int n, int b, int l) {
    return 1;
}

int bkw_lf1(math_t * res, int n, long q, int b, int l, math_t *** T,
            math_t ** aux) {
    size_t i, idx;

    while (1) {
        /* samples from previous layer */
        if (l == 0) {
            return lwe_oracle(res, n, q);
        } else {
            if (!bkw_lf1(aux[0], n, q, b, l - 1, T, aux + 1)) {
                return 0; /* false */
            }
        }

        /* if first elements already 0, returns it */
        if (zero(aux[0], (l - 1) * b, l * b)) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = aux[0][i].value;
            }
            return 1; /* true */
        }

        /* checks if collision */
        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (aux[0][i].value + q - T[l][idx][i].value) % q;
            }
            return 1; /* true */
        }

        /* negates it */
        for (i = 0; i < n + 1; ++i) {
            aux[0][i].value = (q - aux[0][i].value) % q;
        }

        /* checks if collision */
        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (aux[0][i].value + q - T[l][idx][i].value) % q;
            }
            return 1; /* true */
        }

        /* otherwise, stores it, and repeat */
        T[l][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l][idx][i].value = aux[0][i].value;
        }
    }
    return 0; /* should not happen */
}

void bkw_hypo_testing(double ** S, math_t ** F, int d, int m, long q, double sigma,
                      math_t ** aux) {
    int i;
    size_t j, idx;
    double p;

    /* precomputes each value of log(j) */
    if (log_j == NULL) {
        log_j = malloc(q * sizeof(double));
        for (i = 0; i < q; ++i) {
            p = rounded_gaussian_pdf(i, sigma, q, 10*d);
            log_j[i] = log(p) - log((pow(q, d-1) - p)/(double)(pow(q, d) - 1.0));
        }
    }

    /* aux[0][i] = first coordinate of a_i */
    for (i = 0; i < m; ++i) {
        aux[0][i].value = F[i][0].value;
    }

    idx = 1;
    /* foreach v (aux[1]) in Z_q^d, starting at v = (1 0 0 0...) */
    for (aux[1][0].value = 1; aux[1][d - 2].value != q; ) {
        /* foreach sample, update S[v] with a_i - c_i */
        for (i = 0; i < m; ++i) {
            /* S[v]_i = W(aux[0][i] - c_i) */
            S[i][idx] = log_j[(aux[0][i].value - F[i][d - 1].value + q) % q];
            /* aux[0][i] += first coordinate of a_i */
            aux[0][i].value = (aux[0][i].value + F[i][0].value) % q;
        }

        /* increases v by (1 0 0 0...) */
        aux[1][0].value += 1;
        idx++;

        /* when v becomes (q 0 0 0...), increases the next coordinates,
        repeats the process if next coordinates also hit q (max). */
        for (i = 0; aux[1][i].value >= q && i < d - 2 ; ++i) {
            aux[1][i].value = 0;
            /* if there is one coordinate after this one */
            if (i + 1 < d - 1) {
                aux[1][i + 1].value += 1;
                for (j = 0; j < m; ++j) {
                    aux[0][j].value = aux[0][j].value + F[j][i + 1].value;
                }
            } else {
                aux[1][i].value = q; /* to get outside of the loop */
            }
        }
    }
}

void bkw_free_log() {
    free(log_j);
    log_j = NULL;
}
