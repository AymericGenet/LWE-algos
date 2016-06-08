/*
 * bkw_mem_test.c
 *
 *  Created on: May 29, 2016
 *      Author: Aymeric Genet
 */

#include "bkw_mem_test.h"
#include "minunit.h"
#include "../src/bkw_mem.h"
#include "../src/lwe.h"
#include "../src/math.h"
#include "../src/misc.h"
#include <stdio.h>
#include <math.h>

char * test_bkw_mem_lf1() {
    size_t i, j;
    int n, a, d, m;
    long q;
    vec_t res;
    math_t ** aux;
    lwe_t lwe;
    bkw_mem_t bkw;

    /* init  */
    m = 1;
    n = 9, q = 5, d = 3, a = 3;

    secret = malloc(n * sizeof(math_t));
    secret[0] = 1;
    secret[1] = 2;
    secret[2] = 3;
    secret[3] = 4;
    secret[4] = 0;
    secret[5] = 4;
    secret[6] = 3;
    secret[7] = 2;
    secret[8] = 1;
    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    lwe_create(&lwe, n, q, rounded_gaussian, sigma);
    bkw_mem_create(&bkw, lwe, a, d, m);

    res = calloc((n + 1), sizeof(math_t));
    aux = malloc(3*a * sizeof(math_t *));

    for (i = 0; i < a; ++i) {
        for (j = 0; j < 3; ++j) {
            aux[3*i + j] = calloc((n + 1), sizeof(math_t));
        }
    }

    /* draws a sample when l = 0 */
    bkw_mem_lf1(res, bkw, 0, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 1 */
    bkw_mem_lf1(res, bkw, 1, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 2 */
    bkw_mem_lf1(res, bkw, 2, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    for (i = 0; i < a; ++i) {
        printf("\n\t### Table T[%i] state ###\n\n", i);
        rewind(bkw.tables[i]);
        while (read_sample(res, bkw.tables[i], q, n + 1)) {
            printf("\t[ ");
            for (j = 0; j < n + 1; ++j) {
                printf("%lu ", res[j]);
            }
            printf("]\n");
        }
    }
    printf("\n");

    /* frees memory */
    for (i = 0; i < a; ++i) {
        for (j = 0; j < 3; ++j) {
            free(aux[3*i + j]);
        }
    }

    bkw_mem_free(&bkw);

    free(aux);
    free(res);
    free(secret);

    return NULL;
}

char * test_bkw_mem_lf2() {
    size_t i, j;
    int n, a, d, m;
    long q;
    vec_t res;
    math_t ** aux;
    lwe_t lwe;
    bkw_mem_t bkw;

    /* init  */
    m = 1;
    n = 9, q = 5, d = 3, a = 3;

    secret = malloc(n * sizeof(long));
    secret[0] = 1;
    secret[1] = 2;
    secret[2] = 3;
    secret[3] = 4;
    secret[4] = 0;
    secret[5] = 4;
    secret[6] = 3;
    secret[7] = 2;
    secret[8] = 1;
    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    lwe_create(&lwe, n, q, rounded_gaussian, sigma);
    bkw_mem_create(&bkw, lwe, a, d, m);

    res = calloc(n + 1, sizeof(math_t));
    aux = malloc(a * sizeof(math_t *));

    for (i = 0; i < a; ++i) {
        aux[i] = calloc(n + 1, sizeof(math_t));
    }

    /* draws a sample when l = 0 */
    bkw_mem_lf2(res, bkw, 0, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 1 */
    bkw_mem_lf2(res, bkw, 1, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 2 */
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);
    bkw_mem_lf2(res, bkw, 2, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    for (i = 0; i < a; ++i) {
        printf("\n\t### Table T[%i] state ###\n\n", i);
        rewind(bkw.tables[i]);
        while (read_sample(res, bkw.tables[i], q, n + 1)) {
            printf("\t[ ");
            for (j = 0; j < n + 1; ++j) {
                printf("%lu ", res[j]);
            }
            printf("]\n");
        }
    }
    printf("\n");

    /* frees memory */
    for (i = 0; i < a; ++i) {
        free(aux[i]);
    }

    bkw_mem_free(&bkw);

    free(aux);
    free(res);
    free(secret);

    return NULL;
}

char * bkw_mem_all_tests() {
    printf("\n============== BKW Algorithm tests =============\n\n");
    mu_run_test(test_bkw_mem_lf1, "bkw_mem_lf1()");
    mu_run_test(test_bkw_mem_lf2, "bkw_mem_lf2()");
    printf("\n");

    return NULL;
}
