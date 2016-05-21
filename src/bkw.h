/*
 * bkw.h
 *  Implements the BKW algorithm with the following variants for the three
 *  different stages :
 *      1) Sample Reduction   : LF1, LF2, Lazy Modulus Switching
 *      2) Hypothesis Testing : Log-likelihood, FFT
 *      3) Back substitution  : straightforward
 *
 *  Created on: Apr 28, 2016
 *      Author: Aymeric Genet
 */

#include "math.h"

/*
 * Runs the BKW algorithm (to be completed...)
 */

int bkw_algo(math_t * res, int n, int b, int l);

/*
 * Draws a pair (a, c) from the lwe_oracle() such that the components between
 * [(l - 1)*b, l*b] are all zero.
 *
 * Returns 1 (true) if LF1 could sample from the oracle, 0 (false) otherwise. If
 * the function returned 0 (false), the sample in res should be discarded.
 *
 *  Minimum size for auxiliary variable : [l] x [n + 1]
 */

int bkw_lf1(math_t * res, int n, long q, int b, int l, math_t *** T,
            math_t ** aux);

/*
 * Computes the score of every possible vector v in Z_q^b.
 * 
 * The function starts by precomputing the logarithm of every possible noise
 * value. It then proceeds by iterating over all v, starting with
 * v = (1 0 0 0 ...), and computes the value of the scalar product by adding
 * components one by one. Finally, it puts the score of every corresponding
 * vector in the table S.
 *
 *  Source : Albrecht et al, On the complexity of the BKW algorithm on LWE, 2012
 *
 *  Minimum size for auxiliary variable : [2] x [...]
 *    - aux[0] : [m]
 *    - aux[1] : [d - 1]
 */

void bkw_hypo_testing(double ** S, math_t ** F, int d, int m, long q, double sigma,
                      math_t ** aux);

/*
 * Frees the precomputation of the logarithms in bkw_hypo_testing(). This should
 * be called at the end of an execution and if you care about memory leaks...
 */

void bkw_free_log();
