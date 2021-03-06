/*
 * bkw_test.c
 *
 *  Created on: Apr 28, 2016
 *      Author: Aymeric Genet
 */

#include "bkw_test.h"
#include "minunit.h"
#include "../src/bkw.h"
#include "../src/lwe.h"
#include "../src/math.h"
#include "../src/misc.h"

#include <stdio.h>
#include <math.h>

char * test_bkw_lf1() {
    size_t i, j, k;
    int n, a, d, b, m;
    long q, depth;
    double sigma;
    vec_t res;
    math_t ** aux;
    bkw_t * bkw;
    lwe_t * lwe;
    vec_t secret;

    /* init */
    m = 1;
    n = 6, q = 37, d = 1, b = 2, a = 4;
    depth = pow(q, b);

    secret = malloc(n * sizeof(long));
    secret[0] = 1;
    secret[1] = 2;
    secret[2] = 3;
    secret[3] = 4;
    secret[4] = 0;
    secret[5] = 4;

    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    res = calloc((n + 1), sizeof(math_t));
    aux = malloc(a * sizeof(vec_t));

    lwe = malloc(sizeof(lwe_t));
    bkw = malloc(sizeof(bkw_t));
    lwe_create(lwe, n, q, rounded_gaussian, sigma, secret);
    bkw_create(bkw, lwe, a, b, d, m);

    for (i = 0; i < a; ++i) {
        aux[i] = calloc((n + 1), sizeof(math_t));
    }

    /* draws a sample when l = 0 */
    bkw_lf1(res, bkw, 0, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 1 */
    bkw_lf1(res, bkw, 1, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = 2 */
    bkw_lf1(res, bkw, 2, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n");

    /* draws a sample when l = a */
    bkw_lf1(res, bkw, a, aux);

    for (i = 0; i < n + 1; ++i) {
        printf("\tres[%i] = %lu\n", i, res[i]);
    }

    printf("\n\t### Table T[a][q^b] state ###\n\n");
    for (i = 0; i < a; ++i) {
        for (j = 0; j < depth; ++j) {
            if (bkw->tab->first->T[i][j] != NULL) {
                printf("\t[ ");
                for (k = 0; k < n + 1; ++k) {
                    printf("%lu ", bkw->tab->first->T[i][j][k]);
                }
                printf("]\n");
            }
        }
    }
    printf("\n");

    /* frees memory */
    bkw_free(bkw);
    lwe_free(lwe);
    free(bkw);
    free(lwe);
    for (i = 0; i < a; ++i) {
        free(aux[i]);
    }

    free(aux);
    free(res);
    free(secret);

    return NULL;
}

char * test_bkw_lf2() {
    size_t i, j, k, l;
    int n, a, d, b, m;
    long q, depth;
    vec_t res;
    math_t ** aux;
    node_t * curr;
    vec_t ** T;
    bkw_t * bkw;
    lwe_t * lwe;
    vec_t secret;
    double sigma;

    /* init  */
    m = 1;
    n = 6, q = 37, a = 5, b = 1, d = n - (a - 1)*b;
    depth = pow(q, b);

    printf("n = %i, a = %i, b = %i, d = %i\n", n, a, b, d);

    secret = malloc(n * sizeof(long));
    secret[0] = 1;
    secret[1] = 2;
    secret[2] = 3;
    secret[3] = 4;
    secret[4] = 0;
    secret[5] = 4;
    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    lwe = malloc(sizeof(lwe_t));
    bkw = malloc(sizeof(bkw_t));
    lwe_create(lwe, n, q, rounded_gaussian, sigma, secret);
    bkw_create(bkw, lwe, a, b, d, m);

    res = calloc((n + 1), sizeof(math_t));
    aux = malloc(a * sizeof(math_t *));

    for (i = 0; i < a; ++i) {
        aux[i] = calloc((n + 1), sizeof(math_t));
    }

    mem_used = 0;
    lwe_oracle_calls = 0;

    /* draws a sample when l = a */
    for (j = 0; j < 20; ++j) {
        do {
            bkw_lf2(res, bkw, a - 1, aux);
        } while(zero(res, 0, n));
        printf("[ ");
        for (i = 0; i < n + 1; ++i) {
            printf("%lu ", res[i]);
        }
        printf("]\n");
    }

    curr = bkw->tab->first;
    l = 0;
    do {
        T = curr->T;
        printf("\n\t### Linked-Layer %i ###\n\n", l);
        for (i = 0; i < a; ++i) {
            for (j = 0; j < depth; ++j) {
                if (T[i][j] != NULL) {
                    printf("\t[ ");
                    for (k = 0; k < n + 1; ++k) {
                        printf("%lu ", T[i][j][k]);
                    }
                    printf("]\n");
                }
            }
        }
        printf("\n");
        curr = curr->next;
        l++;
    } while (curr != NULL);

    printf("Amount of memory used  : %lu\n", mem_used);
    printf("Amount of oracle calls : %lu\n", lwe_oracle_calls);

    /* frees memory */
    bkw_free(bkw);
    lwe_free(lwe);
    free(bkw);
    free(lwe);

    for (i = 0; i < a; ++i) {
        free(aux[i]);
    }

    free(aux);
    free(res);
    free(secret);

    return NULL;
}

char * test_bkw_distrib() {
    size_t i, j;
    int n, a, b, d, m;
    long noise, q;
    vec_t * res;
    math_t ** aux;
    bkw_t * bkw;
    lwe_t * lwe;
    int * stats;
    double dist;
    vec_t secret;
    double sigma;

    /* init */
    m = 1000;
    n = 6, q = 37, a = 5, d = 2, b = 1;
    stats = calloc(q, sizeof(int));

    secret = malloc(n * sizeof(long));
    for (i = 0; i < n; ++i) {
        read_drandom(&dist);
        secret[i] = dist * (q - 1);
    }
    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    res = malloc(m * sizeof(math_t));
    aux = malloc(a * sizeof(vec_t));

    lwe = malloc(sizeof(lwe_t));
    bkw = malloc(sizeof(bkw_t));
    lwe_create(lwe, n, q, rounded_gaussian, sigma, secret);
    bkw_create(bkw, lwe, a, b, d, m);

    for (i = 0; i < a; ++i) {
        aux[i] = calloc((n + 1), sizeof(math_t));
    }

    for (i = 0; i < m; ++i) {
        res[i] = calloc(n + 1, sizeof(math_t));
    }

    /* collects samples */
    for (i = 0; i < m; ++i) {
        do {
            bkw_lf2(res[i], bkw, a - 1, aux);
        } while(zero(res[i], 0, n));
        noise = res[i][n];
        for (j = 0; j < n; ++j) {
            noise = (q*q + noise - res[i][j]*secret[j]) % q;
        }
        stats[noise]++;
    }

    /* prints distribution */
    printf("\n");
    for (i = 0; i < q; ++i) {
        printf("%i : \t", i);
        for (j = 0; j < 100.0*stats[i]/((double)m); ++j) {
            printf("=");
        }
        dist = 0.0;
        printf(" (%f vs %f)\n", stats[i]/((double) m), rounded_gaussian_pdf(i, sqrt(pow(2, a-2)) * sigma, 1, q));
    }
    printf("\n");

    /* frees memory */
    bkw_free(bkw);
    lwe_free(lwe);
    free(bkw);
    free(lwe);
    for (i = 0; i < a; ++i) {
        free(aux[i]);
    }

    for (i = 0; i < m; ++i) {
        free(res[i]);
    }

    free(stats);
    free(aux);
    free(res);
    free(secret);

    return NULL;
}

char * test_bkw_hypo_testing() {
    size_t i;
    int d, m, a, b, n;
    long q;
    math_t ** aux;
    vec_t * F;
    vec_t v;
    bkw_t * bkw;
    lwe_t * lwe;
    vec_t secret;
    double sigma;

    /* init */
    n = 1, a = 1, b = 1;
    m = 5, q = 5, d = 3, sigma = 1.0;

    secret = calloc(n, sizeof(math_t));
    lwe = malloc(sizeof(lwe_t));
    bkw = malloc(sizeof(bkw_t));
    lwe_create(lwe, n, q, rounded_gaussian, sigma, secret);
    bkw_create(bkw, lwe, a, b, d, m);

    aux = malloc(2 * sizeof(math_t *));
    F = malloc(m * sizeof(vec_t));
    v = calloc(d, sizeof(math_t));

    for (i = 0; i < m; ++i) {
        F[i] = malloc((d + 1) * sizeof(math_t));
    }

    aux[0] = malloc(m * sizeof(math_t));
    aux[1] = calloc(d, sizeof(math_t));

    /* samples */
    F[0][0] = 1; F[0][1] = 2; F[0][2] = 3;
    F[0][3] = (1*1 + 2*2 + 3*3) % q; /* c_0 */

    F[1][0] = 4; F[1][1] = 3; F[1][2] = 2;
    F[1][3] = (4*1 + 3*2 + 2*3) % q; /* c_1 */

    F[2][0] = 2; F[2][1] = 3; F[2][2] = 4;
    F[2][3] = (2*1 + 3*2 + 4*3) % q; /* c_2 */

    F[3][0] = 1; F[3][1] = 4; F[3][2] = 2;
    F[3][3] = (1*1 + 4*2 + 2*3) % q; /* c_3 */

    F[4][0] = 3; F[4][1] = 3; F[4][2] = 4;
    F[4][3] = (3*1 + 3*2 + 4*3) % q; /* c_4 */

    /* runs hypothesis testing */
    bkw_hypo_testing(v, F, bkw, aux);

    /* checks result */
    printf("\tAccording to L-L  v = [ ");
    for (i = 0; i < d; ++i) {
        printf("%lu ", v[i]);
    }
    printf("]\n");
    printf("\tCorrect vector    s = [ 1 2 3 ]\n\n");

    /* frees memory */
    bkw_free(bkw);
    lwe_free(lwe);
    free(bkw);
    free(lwe);
    bkw_free_log();

    for (i = 0; i < m; ++i) {
        free(F[i]);
    }

    free(aux[0]);
    free(aux[1]);

    free(aux);
    free(v);
    free(F);
    free(secret);

    return NULL;
}

char * test_bkw_fft() {
    size_t i;
    int d, m, a, b, n;
    long q;
    double sigma;
    vec_t * F;
    vec_t v;
    bkw_t * bkw;
    lwe_t * lwe;
    vec_t secret;

    /* init */
    n = 1, sigma = 1.0, a = 1, b = 1;
    m = 5, q = 5, d = 3;

    secret = calloc(n, sizeof(math_t));
    lwe = malloc(sizeof(lwe_t));
    bkw = malloc(sizeof(bkw_t));
    lwe_create(lwe, n, q, rounded_gaussian, sigma, secret);
    bkw_create(bkw, lwe, a, b, d, m);

    F = malloc(m * sizeof(vec_t));
    v = calloc(d, sizeof(math_t));

    for (i = 0; i < m; ++i) {
        F[i] = malloc((d + 1) * sizeof(math_t));
    }

    /* samples */
    F[0][0] = 1; F[0][1] = 2; F[0][2] = 3;
    F[0][3] = (1*1 + 2*2 + 3*3) % q; /* c_0 */

    F[1][0] = 4; F[1][1] = 3; F[1][2] = 2;
    F[1][3] = (4*1 + 3*2 + 2*3) % q; /* c_1 */

    F[2][0] = 2; F[2][1] = 3; F[2][2] = 4;
    F[2][3] = (2*1 + 3*2 + 4*3) % q; /* c_2 */

    F[3][0] = 1; F[3][1] = 4; F[3][2] = 2;
    F[3][3] = (1*1 + 4*2 + 2*3) % q; /* c_3 */

    F[4][0] = 3; F[4][1] = 3; F[4][2] = 4;
    F[4][3] = (3*1 + 3*2 + 4*3) % q; /* c_4 */

    /* runs hypothesis testing with fft */
    bkw_fft(v, F, bkw);

    /* checks result */
    printf("\tAccording to FFT  v = [ ");
    for (i = 0; i < d; ++i) {
        printf("%lu ", v[i]);
    }
    printf("]\n");
    printf("\tCorrect vector    s = [ 1 2 3 ]\n\n");

    /* frees memory */
    bkw_free(bkw);
    lwe_free(lwe);
    free(bkw);
    free(lwe);

    for (i = 0; i < m; ++i) {
        free(F[i]);
    }

    free(F);
    free(v);
    free(secret);

    return NULL;
}

char * bkw_all_tests() {
    printf("\n============== BKW Algorithm tests =============\n\n");
    mu_run_test(test_bkw_lf1, "bkw_lf1()");
    mu_run_test(test_bkw_lf2, "bkw_lf2()");
    mu_run_test(test_bkw_distrib, "bkw_distrib()");
    mu_run_test(test_bkw_hypo_testing, "bkw_hypo_testing()");
    mu_run_test(test_bkw_fft, "bkw_fft()");
    printf("\n");

    return NULL;
}
