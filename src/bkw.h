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
#include "fftw3.h"

/* Linked-list node of tables T[a][q^b - 1] */
typedef struct node_t node_t;

struct node_t {
    node_t * next;
    math_t *** T;
};

/*
 * Table structure which allows to keep track of the samples in LF2.
 */

typedef struct {
    node_t * first;
    int * states;
    math_t ** sample;

    int n, b, a, d;
    long q;
} table_t;

void bkw_create_table(table_t * tab, int n, long q, int b, int d);

void bkw_free_table(table_t * tab);

void bkw_create_node(node_t * node, int a, long q, int b);

void bkw_free_node(node_t * node, int a, int q, int b);

/*
 * Runs the BKW algorithm (to be completed...)
 */

int bkw_algo(math_t * res, int n, int b, int l);

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

int bkw_lf1(math_t * res, int n, long q, int b, int l, math_t *** T,
            math_t ** aux);

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

int bkw_lf2(math_t * res, int n, long q, int b, int l, table_t * tab,
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
 * Solves an LWE instance with multi-dimensional fast Fourier transforms.
 */

void bkw_fft(math_t * v, math_t ** F, int d, int m, long q);

/*
 * Frees the precomputation of the logarithms in bkw_hypo_testing(). This should
 * be called at the end of an execution and if you care about memory leaks...
 */

void bkw_free_log();
