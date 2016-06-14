/*
 * lwe_test.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "lwe.h"
#include "misc.h"


void lwe_create(lwe_t * lwe, int n, long q, distribution_t distrib, double sig,
                vec_t secret) {
    size_t i;

    lwe->n = n;
    lwe->q = q;
    lwe->distrib = distrib;
    lwe->sig = sig;
    lwe->secret = malloc(n * sizeof(math_t));
    for (i = 0; i < n; ++i) {
        lwe->secret[i] = secret[i] % q;
    }
}

void lwe_free(lwe_t * lwe) {
    free(lwe->secret);
}

int lwe_oracle(vec_t res, lwe_t * lwe) {
    size_t i;
    double u = 0.0;
    int n = lwe->n;
    long q = lwe->q;

    res[n] = random_sample(lwe->distrib, lwe->sig, q);
    for (i = 0; i < n; ++i) {
        read_drandom(&u);

        res[i] = u * (q - 1);
        res[n] = (res[n] + res[i]*lwe->secret[i]) % q;
    }

    return 1; /* true */
}
