/*
 * lwe_test.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "lwe_test.h"
#include "minunit.h"
#include "../src/lwe.h"
#include "../src/math.h"
#include <stdio.h>
#include <math.h>


char * test_lwe_oracle_predef() {
    int n;
    vec_t pair;
    long q;
    long * s;
    double sig;
    size_t i;

    n = 3;
    q = 97;
    pair = malloc((n + 1) * sizeof(math_t));
    s = malloc(n * sizeof(long));
    s[0] = 35;
    s[1] = 39;
    s[2] = 94;
    sig = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    lwe_oracle_predef(pair, s, n, q, sig);
    printf("\n");
    for (i = 0; i < n; ++i) {
        printf("\ta[%i] = %lu\n", i, pair[i]);
    }
    printf("\tc = %lu\n", pair[n]);

    free(pair);
    free(s);

    return NULL;
}

char * test_lwe_oracle() {
    int n;
    vec_t pair;
    long q;
    size_t i;

    n = 3;
    q = 97;
    pair = malloc((n + 1) * sizeof(math_t));
    secret = malloc(n * sizeof(long));
    secret[0] = 35;
    secret[1] = 39;
    secret[2] = 94;
    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    lwe_oracle_predef(pair, secret, n, q, sigma);
    printf("\n");
    for (i = 0; i < n; ++i) {
        printf("\ta[%i] = %lu\n", i, pair[i]);
    }
    printf("\tc = %lu\n", pair[n]);

    free(pair);
    free(secret);

    return NULL;
}

char * lwe_all_tests() {
    printf("\n=============== LWE Oracle tests ===============\n\n");
    mu_run_test(test_lwe_oracle_predef, "lwe_oracle_predef()");
    mu_run_test(test_lwe_oracle, "lwe_oracle()");
    printf("\n");

    return NULL;
}
