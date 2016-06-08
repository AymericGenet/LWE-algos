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
} lwe_t;

extern vec_t secret;
extern double sigma;


void lwe_create(lwe_t * lwe, int n, long q, distribution_t distrib, double sig);

/*
 * Takes secret s and puts in res an array containing the pair (a, c) such that
 * a is uniformly in Z_q^n, and c = <a, s> + e where e follows a distribution D.
 *
 * Returns 1 (true) if query was successful, 0 (false) otherwise.
*/

int lwe_oracle_predef(vec_t res, vec_t s, lwe_t lwe);

/*
 * Draws a pair (a, c) such that a is uniformly in Z_q^n, and c = <a, s> + e
 * where e follows a distribution D.
 *
 * The function defines an oracle of a real case scenario, it can be linked to
 * a file.
 * 
 * Returns 1 (true) if query was successful, 0 (false) otherwise.
 */

int lwe_oracle(vec_t res, lwe_t lwe);

#endif /* LWE_H_ */
