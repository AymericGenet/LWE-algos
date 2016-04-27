/*
 * lwe_oracle_test.c : Unit tests for LWE oracle.
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "lwe_oracle_test.h"
#include "minunit.h"
#include "../src/lwe_oracle.h"
#include "../src/math.h"
#include <stdio.h>
#include <math.h>

char * test_lwe_oracle() {
    int n;
    math_t * pair;
    long q;
    long * s;
    double sigma;
    size_t i;

    n = 3;
    q = 97;
    pair = malloc((n + 1) * sizeof(math_t));
    s = malloc(n * sizeof(long));
    s[0] = 35;
    s[1] = 39;
    s[2] = 94;
    sigma = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    lwe_oracle(pair, s, n, q, sigma);
    printf("\n");
    for (i = 0; i < n; ++i) {
        printf("\ta[%i] = %lu\n", i, pair[i]);
    }
    printf("\tc = %lu\n", pair[n]);

    return NULL;
}

char * lwe_oracle_all_tests() {
    printf("\n=============== LWE Oracle tests ===============\n\n");
    mu_run_test(test_lwe_oracle, "lwe_oracle()");
    printf("\n");

    return NULL;
}
