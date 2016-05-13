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

char * test_bkw_hypo_testing() {
    size_t i, j;
    int d, m;
    long q, depth;
    math_t ** aux;
    math_t ** S;
    math_t ** F;

    m = 2, q = 5, d = 4;
    depth = 125; /* q^(d - 1) */

    S = malloc(m * sizeof(math_t *));
    aux = malloc(m * sizeof(math_t *));
    F = malloc(m * sizeof(math_t *));

    for (i = 0; i < m; ++i) {
        S[i] = malloc(depth * sizeof(math_t));
        aux[i] = malloc(d * sizeof(math_t));
        F[i] = malloc(d * sizeof(math_t));
        for (j = 0; j < d; ++j) {
            aux[i][j].value = 0;
        }
    }

    F[0][0].value = 1; F[0][1].value = 2; F[0][2].value = 3; F[0][3].value = 4;
    F[1][0].value = 4; F[1][1].value = 3; F[1][2].value = 2; F[1][3].value = 1;
/*    F[2][0].value = 9; F[2][1].value = 10; F[2][2].value = 11; F[2][3].value = 12;
    F[3][0].value = 13; F[3][1].value = 14; F[3][2].value = 15; F[3][3].value = 16;
    F[4][0].value = 17; F[4][1].value = 18; F[4][2].value = 19; F[4][3].value = 20;*/

    bkw_hypo_testing(S, F, d, m, q, aux);

    printf("Done ! \n");

    for (i = 0; i < depth; ++i) {
        printf("[ ");
        for (j = 0; j < m; ++j) {
            printf("%lu ", S[j][i].value);
        }
        printf("]\n");
    }

    return NULL;
}

char * bkw_all_tests() {
    printf("\n============== BKW Algorithm tests =============\n\n");
    mu_run_test(test_bkw_lf1, "bkw_lf1()");
    mu_run_test(test_bkw_hypo_testing, "bkw_hypo_testing()");
    printf("\n");

    return NULL;
}
