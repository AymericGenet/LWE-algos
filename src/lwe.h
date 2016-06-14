/*
 * lwe_test.h
 *  Provides LWE oracle functions that will be used in the BKW algorithm.
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#ifndef LWE_H_
#define LWE_H_

#define NOISE_FUNCTION(x, q) rounded_gaussian(x, q)

#include "math.h"

typedef struct {
    int n;
    long q;
    distribution_t distrib;
    double sig;
    vec_t secret;
} lwe_t;

extern vec_t secret;
extern double sigma;


void lwe_create(lwe_t * lwe, int n, long q, distribution_t distrib, double sig,
                vec_t secret);

void lwe_free(lwe_t * lwe);

/*
 * Draws a pair (a, c) according to LWE instance.
 * 
 * Returns 1 (true) if query was successful, 0 (false) otherwise.
 */

int lwe_oracle(vec_t res, lwe_t * lwe);

#endif /* LWE_H_ */
