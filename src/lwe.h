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

extern long * secret;
extern double sigma;

/*
 * Takes secret s and puts in res an array containing the pair (a, c) such that
 * a is uniformly in Z_q^n, and c = <a, s> + e where e follows a distribution D.
 *
 * Returns 1 (true) if query was successful, 0 (false) otherwise.
*/

int lwe_oracle_predef(math_t * res, long * s, int n, long q, double sig);

/*
 * Draws a pair (a, c) such that a is uniformly in Z_q^n, and c = <a, s> + e
 * where e follows a distribution D.
 *
 * The function defines an oracle of a real case scenario, it can be linked to
 * a file.
 * 
 * Returns 1 (true) if query was successful, 0 (false) otherwise.
 */

int lwe_oracle(math_t * res, int n, long q);

#endif /* LWE_H_ */