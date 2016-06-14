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


char * test_lwe_oracle() {
    int n;
    vec_t pair;
    long q;
    vec_t s;
    double sig;
    size_t i;
    lwe_t * lwe;

    n = 3;
    q = 97;
    pair = malloc((n + 1) * sizeof(math_t));
    s = malloc(n * sizeof(math_t));
    s[0] = 35;
    s[1] = 39;
    s[2] = 94;
    sig = q/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    lwe = malloc(sizeof(lwe_t));
    lwe_create(lwe, n, q, rounded_gaussian, sig, s);

    lwe_oracle(pair, lwe);
    printf("\n");
    for (i = 0; i < n; ++i) {
        printf("\ta[%i] = %lu\n", i, pair[i]);
    }
    printf("\tc = %lu\n", pair[n]);

    lwe_free(lwe);
    free(lwe);
    free(pair);
    free(s);

    return NULL;
}

char * lwe_all_tests() {
    printf("\n=============== LWE Oracle tests ===============\n\n");
    mu_run_test(test_lwe_oracle, "lwe_oracle()");
    printf("\n");

    return NULL;
}
