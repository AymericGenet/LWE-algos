/*
 * lwe_test.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "lwe.h"
#include "math.h"

long * secret;
double sigma;

int lwe_oracle_predef(vec_t res, long * s, int n, long q, double sig) {
    size_t i;

    res[n] = rounded_gaussian(sig, q);
    for (i = 0; i < n; ++i) {
        read_random(res + i);

        res[i] %= q;
        res[n] = (res[n] + res[i] * s[i]) % q;
    }

    return 1; /* true */
}

int lwe_oracle(vec_t res, int n, long q) {
    return lwe_oracle_predef(res, secret, n, q, sigma);
}
