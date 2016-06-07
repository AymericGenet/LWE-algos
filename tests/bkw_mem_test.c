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
    int n, b, a;
    long q;
    vec_t res;
    math_t ** aux;
    FILE ** tables;
    char path[128];

    /* init  */
    n = 9, q = 5, b = 3;
    a = n/b;

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

    res = calloc((n + 1), sizeof(math_t));
    aux = malloc(3*a * sizeof(math_t *));
    tables = malloc(a * sizeof(FILE *));

    for (i = 0; i < a; ++i) {
        for (j = 0; j < 3; ++j) {
            aux[3*i + j] = calloc((n + 1), sizeof(math_t));
        }
        sprintf(path, "tables/lf1.T.%i.txt", i);
        open_table(&(tables[i]), path);
    }

    /* draws a sample when l = 0 */
    bkw_mem_lf1(res, n, q, b, 0, tables, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 1 */
    bkw_mem_lf1(res, n, q, b, 1, tables, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 2 */
    bkw_mem_lf1(res, n, q, b, 2, tables, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    for (i = 0; i < a; ++i) {
        printf("\n\t### Table T[%i] state ###\n\n", i);
        rewind(tables[i]);
        while (read_sample(res, tables[i], q, n + 1)) {
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
        close_table(&(tables[i]));
    }

    free(aux);
    free(tables);
    free(res);
    free(secret);

    return NULL;
}

char * test_bkw_mem_lf2() {
    size_t i, j;
    bkw_mem_t * bkw;
    vec_t res;
    math_t ** aux;
    char path[128];

    /* init  */
    bkw = malloc(sizeof(bkw_mem_t));
    bkw->n = 9;
    bkw->q = 5;
    bkw->b = 3;
    bkw->d = 2;
    bkw->a = bkw->n/bkw->b;
    bkw->tables = malloc(bkw->a * sizeof(FILE *));
    bkw->sample_pos = calloc(bkw->a, sizeof(vec_t));
    bkw->sample_neg = calloc(bkw->a, sizeof(vec_t));

    secret = malloc(bkw->n * sizeof(long));
    secret[0] = 1;
    secret[1] = 2;
    secret[2] = 3;
    secret[3] = 4;
    secret[4] = 0;
    secret[5] = 4;
    secret[6] = 3;
    secret[7] = 2;
    secret[8] = 1;
    sigma = bkw->q/(sqrt(2 * PI_VAL * bkw->n) * log(bkw->n) * log(bkw->n));

    res = calloc(bkw->n + 1, sizeof(math_t));
    aux = malloc(bkw->a * sizeof(math_t *));

    for (i = 0; i < bkw->a; ++i) {
        aux[i] = calloc(bkw->n + 1, sizeof(math_t));
        sprintf(path, "tables/lf2.T.%i.txt", i);
        open_table(&(bkw->tables[i]), path);
    }

    /* draws a sample when l = 0 */
    bkw_mem_lf2(res, 0, bkw, aux);

    for (i = 0; i < bkw->n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 1 */
    bkw_mem_lf2(res, 1, bkw, aux);

    for (i = 0; i < bkw->n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 2 */
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);
    bkw_mem_lf2(res, 2, bkw, aux);

    for (i = 0; i < bkw->n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    for (i = 0; i < bkw->a; ++i) {
        printf("\n\t### Table T[%i] state ###\n\n", i);
        rewind(bkw->tables[i]);
        while (read_sample(res, bkw->tables[i], bkw->q, bkw->n + 1)) {
            printf("\t[ ");
            for (j = 0; j < bkw->n + 1; ++j) {
                printf("%lu ", res[j]);
            }
            printf("]\n");
        }
    }
    printf("\n");

    /* frees memory */
    for (i = 0; i < bkw->a; ++i) {
        free(aux[i]);
        close_table(&(bkw->tables[i]));
        free(bkw->sample_pos[i]);
        free(bkw->sample_neg[i]);
    }

    free(bkw->sample_pos);
    free(bkw->sample_neg);
    free(bkw->tables);
    free(bkw);
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
