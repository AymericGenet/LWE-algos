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
