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

#ifndef BKW_H_
#define BKW_H_

#include "lwe.h"
#include "math.h"
#include "fftw3.h"

/* Linked-list node of tables T[a][q^b - 1] */
typedef struct node_t node_t;

extern unsigned long lwe_oracle_calls;
extern unsigned long mem_used;

struct node_t {
    node_t * next;
    vec_t ** T;
};

/*
 * Table structure which allows to keep track of the samples in LF2.
 */

typedef struct {
    node_t * first;
    int * states;
    vec_t * sample;
} table_t;

/*
 * BKW instance structure.
 */

typedef struct {
    lwe_t * lwe;
    int a, b, d;
    long m;
    table_t * tab;
} bkw_t;

void bkw_create(bkw_t * bkw, lwe_t * lwe, int a, int b, int d, long m);

void bkw_free(bkw_t * bkw);

void bkw_create_node(node_t * node, int a, long q, int b);

void bkw_free_node(node_t * node, int a, int q, int b);

/*
 * Draws a pair (a, c) from the lwe_oracle() such that the components between
 * [(l - 1)*b, l*b] are all zero, according to LF1.
 *
 * The LF1 algorithm samples the oracle and stores the results in T[l] until the
 * components from [(l - 1)*b, l*b] collides with a previously stored sample. In
 * this case, it outputs the difference between those two. The algorithm also
 * checks for a collision with the negation of the sample, but only stores one
 * of the two.
 *
 * Returns 1 (true) if LF1 could sample from the oracle, 0 (false) otherwise. If
 * the function returned 0 (false), the sample in res should be discarded.
 *
 *  Minimum size for auxiliary variable : [l] x [n + 1]
 */

int bkw_lf1(vec_t res, bkw_t * bkw, int l, math_t ** aux);

/*
 * Draws a pair (a, c) from the lwe_oracle() such that the components between
 * [(l - 1)*b, l*b] are all zero, according to LF2.
 *
 * The LF2 improvement stores every sample, even if a sample collides with one
 * stored in the table T. It makes use of the linked-list structure table_t to
 * append such new sample. 
 *
 * Returns 1 (true) if LF2 could sample from the oracle, 0 (false) otherwise. If
 * the function returned 0 (false), the sample in res should be discarded.
 *
 *  Minimum size for auxiliary variable : [l] x [n + 1]
 */

int bkw_lf2(vec_t res, bkw_t * bkw, int l, math_t ** aux);

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

void bkw_hypo_testing(vec_t v, vec_t * F, bkw_t * bkw, math_t ** aux);

/*
 * Solves an LWE instance with multi-dimensional fast Fourier transforms.
 */

void bkw_fft(vec_t v, vec_t * F, bkw_t * bkw);

/*
 * Frees the precomputation of the logarithms in bkw_hypo_testing(). This should
 * be called at the end of an execution and if you care about memory leaks...
 */

void bkw_free_log();

#endif /* BKW_H_ */
