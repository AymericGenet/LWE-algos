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
    size_t i, j, k;
    int n, b, a;
    long q, depth;
    math_t * res;
    math_t ** aux;
    math_t *** T;

    /* init  */
    n = 9, q = 5, b = 3;
    a = n/b;
    depth = pow(q, b);

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
    aux = malloc(a * sizeof(math_t *));
    T = malloc(a * sizeof(math_t **));

    for (i = 0; i < a; ++i) {
        aux[i] = calloc((n + 1), sizeof(math_t));
        T[i] = calloc(depth, sizeof(math_t *));
    }

    /* draws a sample when l = 0 */
    bkw_lf1(res, n, q, b, 0, T, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i].value);
    }

    printf("\n");

    /* draws a sample when l = 1 */
    bkw_lf1(res, n, q, b, 1, T, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i].value);
    }

    printf("\n");

    /* draws a sample when l = 2 */
    bkw_lf1(res, n, q, b, 2, T, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i].value);
    }

    printf("\n\t### Table T[a][q^b] state ###\n\n");
    for (i = 0; i < a; ++i) {
        for (j = 0; j < depth; ++j) {
            if (T[i][j] != NULL) {
                printf("\t[ ");
                for (k = 0; k < n + 1; ++k) {
                    printf("%lu ", T[i][j][k].value);
                }
                printf("]\n");
            }
        }
    }
    printf("\n");

    /* frees memory */
    for (i = 0; i < a; ++i) {
        free(aux[i]);
        for (j = 0; j < depth; ++j) {
            free(T[i][j]);
        }
        free(T[i]);
    }

    free(aux);
    free(T);
    free(res);
    free(secret);

    return NULL;
}

char * test_bkw_lf2() {
    size_t i, j, k, l;
    int n, b, a, d;
    long q, depth;
    math_t * res;
    math_t ** aux;
    table_t * tab;
    node_t * curr;
    math_t *** T;

    /* init  */
    n = 9, q = 5, b = 3, d = 2;
    a = n/b;
    depth = pow(q, b);

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
    aux = malloc(a * sizeof(math_t *));

    for (i = 0; i < a; ++i) {
        aux[i] = calloc((n + 1), sizeof(math_t));
    }

    tab = malloc(sizeof(table_t));
    bkw_create_table(tab, n, q, b, d);

    /* draws a sample when l = 0 */
    bkw_lf2(res, n, q, b, 0, tab, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i].value);
    }

    printf("\n");

    /* draws a sample when l = 1 */
    bkw_lf2(res, n, q, b, 1, tab, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i].value);
    }

    printf("\n");

    /* draws a sample when l = 2 */
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);
    bkw_lf2(res, n, q, b, 2, tab, aux);

    for (i = 0; i < n; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i].value);
    }

    curr = tab->first;
    l = 0;
    do {
        T = curr->T;
        printf("\n\t### Linked-Layer %i ###\n\n", l);
        for (i = 0; i < a; ++i) {
            for (j = 0; j < depth; ++j) {
                if (T[i][j] != NULL) {
                    printf("\t[ ");
                    for (k = 0; k < n + 1; ++k) {
                        printf("%lu ", T[i][j][k].value);
                    }
                    printf("]\n");
                }
            }
        }
        printf("\n");
        curr = curr->next;
        l++;
    } while (curr != NULL);

    /* frees memory */
    bkw_free_table(tab);

    for (i = 0; i < a; ++i) {
        free(aux[i]);
    }

    free(aux);
    free(res);
    free(secret);

    return NULL;
}

char * test_bkw_hypo_testing() {
    size_t i, j, k;
    int d, m;
    long q, depth;
    math_t ** aux;
    math_t ** F;
    double ** S;
    double best, sum;

    /* init */
    sigma = 1;
    m = 5, q = 5, d = 4;
    depth = 125; /* q^(d - 1) */

    S = malloc(m * sizeof(double *));
    aux = malloc(2 * sizeof(math_t *));
    F = malloc(m * sizeof(math_t *));

    for (i = 0; i < m; ++i) {
        S[i] = calloc(depth, sizeof(double));
        F[i] = malloc(d * sizeof(math_t));
    }

    aux[0] = malloc(m * sizeof(math_t));
    aux[1] = calloc((d - 1), sizeof(math_t));

    /* samples */
    F[0][0].value = 1; F[0][1].value = 2; F[0][2].value = 3;
    F[0][3].value = (1*1 + 2*2 + 3*3) % q; /* c_0 */

    F[1][0].value = 4; F[1][1].value = 3; F[1][2].value = 2;
    F[1][3].value = (4*1 + 3*2 + 2*3) % q; /* c_1 */

    F[2][0].value = 2; F[2][1].value = 3; F[2][2].value = 4;
    F[2][3].value = (2*1 + 3*2 + 4*3) % q; /* c_2 */

    F[3][0].value = 1; F[3][1].value = 4; F[3][2].value = 2;
    F[3][3].value = (1*1 + 4*2 + 2*3) % q; /* c_3 */

    F[4][0].value = 3; F[4][1].value = 3; F[4][2].value = 4;
    F[4][3].value = (3*1 + 3*2 + 4*3) % q; /* c_4 */

    /* runs hypothesis testing */
    bkw_hypo_testing(S, F, d, m, q, sigma, aux);

    /* looks for the best candidate */
    best = 0.0;
    k = 0;
    for (i = 1; i < depth; ++i) {
        printf("\t%i : [ ", i);
        sum = 0.0;
        for (j = 0; j < m; ++j) {
            printf("%.5f ", S[j][i]);
            sum += S[j][i];
        }
        if (sum >= best) {
            k = i;
            best = sum;
        }
        printf("]\n");
    }
    printf("\n");

    printf("\t best score at : %i (%.5f)\n", k, best);

    /* compares with the right secret */
    best = 0.0;
    k = 1 + 2*q + 3*q*q;
    for (j = 0; j < m; ++j) {
        best += S[j][k];
    }

    printf("\t correct was : %i (%.5f)\n\n", k, best);

    /* frees memory */
    bkw_free_log();

    for (i = 0; i < m; ++i) {
        free(S[i]);
        free(F[i]);
    }

    free(aux[0]);
    free(aux[1]);

    free(aux);
    free(S);
    free(F);

    return NULL;
}

char * test_bkw_fft() {
    size_t i;
    int d, m;
    long q;
    math_t ** F;
    math_t * v;

    /* init */
    m = 5, q = 5, d = 4;

    F = malloc(m * sizeof(math_t *));
    v = calloc((d - 1), sizeof(math_t));

    for (i = 0; i < m; ++i) {
        F[i] = malloc(d * sizeof(math_t));
    }

    /* samples */
    F[0][0].value = 1; F[0][1].value = 2; F[0][2].value = 3;
    F[0][3].value = (1*1 + 2*2 + 3*3) % q; /* c_0 */

    F[1][0].value = 4; F[1][1].value = 3; F[1][2].value = 2;
    F[1][3].value = (4*1 + 3*2 + 2*3) % q; /* c_1 */

    F[2][0].value = 2; F[2][1].value = 3; F[2][2].value = 4;
    F[2][3].value = (2*1 + 3*2 + 4*3) % q; /* c_2 */

    F[3][0].value = 1; F[3][1].value = 4; F[3][2].value = 2;
    F[3][3].value = (1*1 + 4*2 + 2*3) % q; /* c_3 */

    F[4][0].value = 3; F[4][1].value = 3; F[4][2].value = 4;
    F[4][3].value = (3*1 + 3*2 + 4*3) % q; /* c_4 */

    /* runs hypothesis testing with fft */
    bkw_fft(v, F, d, m, q);

    /* checks result */
    printf("\tAccording to FFT  v = [ ");
    for (i = 0; i < d - 1; ++i) {
        printf("%lu ", v[i].value);
    }
    printf("]\n");
    printf("\tCorrect vector    s = [ 1 2 3 ]\n\n");

    /* frees memory */
    for (i = 0; i < m; ++i) {
        free(F[i]);
    }

    free(F);
    free(v);

    return NULL;
}

char * bkw_all_tests() {
    printf("\n============== BKW Algorithm tests =============\n\n");
    mu_run_test(test_bkw_lf1, "bkw_lf1()");
    mu_run_test(test_bkw_lf2, "bkw_lf2()");
    mu_run_test(test_bkw_hypo_testing, "bkw_hypo_testing()");
    mu_run_test(test_bkw_fft, "bkw_fft()");
    printf("\n");

    return NULL;
}
