/*
 * bkw_test.c
 *
 *  Created on: Apr 28, 2016
 *      Author: Aymeric Genet
 */

#include "bkw_test.h"
#include "minunit.h"
#include "../src/bkw.h"
#include "../src/lwe_oracle.h"
#include "../src/math.h"
#include <stdio.h>
#include <math.h>

char * test_bkw_lf1() {
    size_t i, j;
    int n, b, a;
    long q, depth;
    math_t * res;
    math_t ** aux;
    math_t *** T;

    n = 9, q = 97, b = 3;
    a = n/b;
    depth = 2*456336; /* (q^b - 1)/2 */

    secret = malloc(n * sizeof(long));
    secret[0] = 1;
    secret[1] = 2;
    secret[2] = 3;
    secret[3] = 4;
    secret[4] = 5;
    secret[5] = 6;
    secret[6] = 7;
    secret[7] = 8;
    secret[8] = 9;
    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    res = malloc((n + 1) * sizeof(math_t));
    aux = malloc(a * sizeof(math_t *));
    T = malloc(a * sizeof(math_t **));

    for (i = 0; i < a; ++i) {
        aux[i] = malloc((n + 1) * sizeof(math_t));
        T[i] = malloc(depth * sizeof(math_t *));

        for (j = 0; j < depth; ++j) {
            T[i][j] = NULL;
        }
    }

    bkw_lf1(res, n, q, b, 0, T, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %i\n", i, res[i].value);
    }

    printf("\n");

    bkw_lf1(res, n, q, b, 1, T, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %i\n", i, res[i].value);
    }

    printf("\n");

    bkw_lf1(res, n, q, b, 2, T, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %i\n", i, res[i].value);
    }

    return NULL;
}

char * bkw_all_tests() {
    printf("\n============== BKW Algorithm tests =============\n\n");
    mu_run_test(test_bkw_lf1, "bkw_lf1()");
    printf("\n");

    return NULL;
}
