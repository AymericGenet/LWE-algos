/*
 * lwe_test.h
 *  Provides LWE oracle functions that will be used in the BKW algorithm.
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#ifndef LWE_H_
#define LWE_H_

#include "math.h"

/*
 * LWE instance data structure.
 */

typedef struct {
    int n;
    long q;
    distribution_t distrib;
    double sig;
    vec_t secret;
} lwe_t;

/*
 * Creates an LWE instance.
 *
 * @param lwe The created LWE instance
 * @param n The size of the secret
 * @param q The modulus of the ring Z_q
 * @param distrib The noise distribution chi
 * @param sig The corresponding standard deviation for the noise distribution
 * @param secret The secret vector
 */

void lwe_create(lwe_t * lwe, int n, long q, distribution_t distrib, double sig,
                vec_t secret);

/*
 * Frees memory allocated by an LWE instance.
 *
 * @param lwe The LWE instance to be deleted
 */

void lwe_free(lwe_t * lwe);

/*
 * Draws a pair (a, c) according to a LWE instance.
 *
 * @param res The resulting oracle sample
 * @param lwe The corresponding LWE instance
 * @return 1 (true) if query was successful, 0 (false) otherwise
 */

int lwe_oracle(vec_t res, lwe_t * lwe);

#endif /* LWE_H_ */
