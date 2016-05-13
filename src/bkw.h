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
 *
 */

int bkw_algo(math_t * res, int n, int b, int l);

/*
 * 
 */

int bkw_lf1(math_t * res, int n, long q, int b, int l, math_t *** T,
            math_t ** aux);

/*
 *
 */

int bkw_hypo_testing(math_t ** S, math_t ** F, int d, int m, long q,
                     math_t ** aux);
