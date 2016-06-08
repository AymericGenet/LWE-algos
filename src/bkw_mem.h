/*
 * bkw_mem.h
 *  Implements the BKW algorithm with the following variants for the three
 *  different stages :
 *      1) Sample Reduction   : LF1, LF2, Lazy Modulus Switching
 *      2) Hypothesis Testing : Log-likelihood, FFT
 *      3) Back substitution  : straightforward
 *
 *  This version stores samples on disk instead on physical memory.
 *
 *  Created on: May 28, 2016
 *      Author: Aymeric Genet
 */

#ifndef BKW_MEM_H_
#define BKW_MEM_H_

#include "bkw.h"
#include "math.h"
#include <stdio.h>

#define BKW_MEM_PATH "tables/T.%i.txt"

typedef struct {
    lwe_t lwe;
    FILE ** tables;
    vec_t * sample_pos;
    vec_t * sample_neg;

    int a, b, d, m;
} bkw_mem_t;


void bkw_mem_create(bkw_mem_t * bkw, lwe_t lwe, int a, int d, int m);

void bkw_mem_free(bkw_mem_t * bkw);

/*
 * Draws a pair (a, c) from the lwe_oracle() such that the components between
 * [(l - 1)*b, l*b] are all zero, according to LF1.
 *
 * The LF1 algorithm samples the oracle and stores the results in "T_l.txt"
 * until the components from [(l - 1)*b, l*b] collides with a previously stored
 * sample. In this case, it outputs the difference between those two. The
 * algorithm also checks for a collision with the negation of the sample, but
 * only stores one of the two.
 *
 * Returns 1 (true) if LF1 could sample from the oracle, 0 (false) otherwise. If
 * the function returned 0 (false), the sample in res should be discarded.
 *
 *  Minimum size for auxiliary variable : [l] x [n + 1]
 */

int bkw_mem_lf1(vec_t res, bkw_mem_t bkw, int l, math_t ** aux);

int bkw_mem_lf2(vec_t res, bkw_mem_t bkw, int l, math_t ** aux);

#endif /* BKW_MEM_H_ */
