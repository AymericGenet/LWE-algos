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


int bkw_algo(math_t * res, int n, int b, int l) {
    return 1;
}

int bkw_lf1(math_t * res, int n, long q, int b, int l, math_t *** T,
            math_t ** aux) {
    size_t i, idx;

    while (1) {
        if (l == 0) {
            return lwe_oracle(res, n, q);
        } else {
            if (!bkw_lf1(aux[0], n, q, b, l - 1, T, aux + 1)) {
                return 0; /* false */
            }
        }

        if (zero(aux[0], (l - 1) * b, l * b)) {
            for (i = 0; i < n; ++i) {
                res[i].value = aux[0][i].value;
            }
            return 1; /* true */
        }

        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n; ++i) {
                res[i].value = (aux[0][i].value - T[l][idx][i].value) % q;
            }

            return 1;
        }

        for (i = 0; i < n; ++i) {
            aux[0][i].value = (q - aux[0][i].value) % q;
        }

        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n; ++i) {
                res[i].value = (aux[0][i].value - T[l][idx][i].value) % q;
            }

            return 1;
        }

        T[l][idx] = malloc(n * sizeof(math_t));
        for (i = 0; i < n; ++i) {
            T[l][idx][i].value = aux[0][i].value;
        }
    }
    return 0; /* should not happen */
}

int bkw_hypo_testing(math_t ** S, math_t ** F, int d, int m, long q,
                     math_t ** aux) {
    int i;
    size_t j, idx;

    for (i = 0; i < m; ++i) {
        aux[0][i].value = F[i][d - 2].value;
    }
    for (aux[1][d - 2].value = 1; aux[1][0].value != q; ) {
        for (i = 0; i < m; ++i) {
            idx = index(aux[1], q, 0, d - 1);
            S[i][idx].value = (aux[0][i].value - F[i][d - 1].value) % q;
            aux[0][i].value = (aux[0][i].value + F[i][d - 2].value) % q;
        }
        aux[1][d - 2].value += 1;
        for (i = d - 2; i >= 0 && aux[1][i].value >= q; --i) {
            aux[1][i].value = 0;
            if (i - 1 >= 0) {
                aux[1][i - 1].value += 1;
                for (j = 0; j < m; ++j) {
                    aux[0][j].value = aux[0][j].value + F[j][i - 1].value;
                }
            } else {
                aux[1][i].value = q;
            }
        }
    }

    return 1; /* true */
}
