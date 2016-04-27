/*
 * lwe_oracle_test.h
 *  Provides LWE oracle functions that will be used in the BKW algorithm.
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#ifndef LWE_ORACLE_H_
#define LWE_ORACLE_H_

#define NOISE_FUNCTION(x) uniform(x)

#include "math.h"

/*
 * Takes secret s and puts in res an array containing the pair (a, c) such that
 * a is uniformly in Z_q^n, and c = <a, s> + e where e follows a distribution D.
 * Returns 1 (true) if query was successful, 0 (false) otherwise.
*/

int lwe_oracle(math_t * res, long * s, int n, long q, double sigma);

#endif /* LWE_ORACLE_H_ */
