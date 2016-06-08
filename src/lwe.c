/*
 * lwe_test.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "lwe.h"
#include "math.h"

vec_t secret;
double sigma;

void lwe_create(lwe_t * lwe, int n, long q, distribution_t distrib, double sig) {
    lwe->n = n;
    lwe->q = q;
    lwe->distrib = distrib;
    lwe->sig = sig;
}

int lwe_oracle_predef(vec_t res, vec_t s, lwe_t lwe) {
    size_t i;
    double u = 0.0;
    int n = lwe.n;
    long q = lwe.q;

    res[n] = random_sample(lwe.distrib, lwe.sig, q);
    for (i = 0; i < n; ++i) {
        read_drandom(&u);

        res[i] = u*q == q ? 0 : u*q; /* because C is really dumb, sometimes */
        res[n] = (res[n] + res[i]*s[i]) % q;
    }

    return 1; /* true */
}

int lwe_oracle(vec_t res, lwe_t lwe) {
    return lwe_oracle_predef(res, secret, lwe);
}
