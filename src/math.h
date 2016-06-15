/*
 * maths.h : Definitions for mathematical operations.
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#ifndef MATH_H_
#define MATH_H_

#include <stdlib.h>

#define PI_VAL 3.14159265358979323846

/* Bounds for wrapping the Gaussian random variable around Z_q */
#define GAUSSIAN_BOUND 10

/* Represents an element of Z_q */
typedef unsigned long math_t;

/* Represents a vector of elements of Z_q */
typedef math_t * vec_t;

/* Different choices of distribution, according to types of LWE noise */
#define DISTRIB_NUMBER 3

typedef enum {
    discrete_gaussian,
    rounded_gaussian,
    uniform
} distribution_t;

/* Access to distributions pdf according to distribution_t */
extern double (*distributions[DISTRIB_NUMBER])(long, double, int, long);

/*
 * Draws a random sample according to the rejection sampling method.
 *
 * @param distrib Distribution according to which it will sample
 * @param sigma The corresponding standard deviation of the distribution
 * @param q The modulus of the ring Z_q
 * @return A sample x drawn at random according to distrib.
 */

long random_sample(distribution_t distrib, double sigma, long q);

/*
 * Computes the probability density function of the discrete Gaussian
 * distribution.
 *
 * @param x The element at which the pdf is evaluated
 * @param sigma The standard deviation of the Gaussian distribution
 * @param q The modulus of the ring Z_q
 * @return The value for the cdf of the rounded Gaussian at x
 */

double discrete_gaussian_pdf(long x, double sigma, int n, long q);

/*
 * Computes the cumulative distribution function of the rounded Gaussian
 * distribution.
 *
 * @param x The element at which the cdf is evaluated
 * @param sigma The standard deviation of the Gaussian distribution
 * @return The value for the cdf of the rounded Gaussian at x
 */

double rounded_gaussian_cdf(double x, double sigma);

/*
 * Computes the probability density function of the rounded Gaussian
 * distribution.
 *
 * @param x The element at which the pdf is evaluated
 * @param sigma The standard deviation of the Gaussian distribution
 * @param q The modulus of the ring Z_q
 * @param k The number of samples for the real Gaussian distribution
 * @return The value for the cdf of the rounded Gaussian at x
 */

double rounded_gaussian_pdf(long x, double sigma, int n, long q);

/*
 * Computes the probability density function of the uniform distribution between
 * [-beta, beta], where
 *      beta = (sqrt(12*sigma*sigma + 1) - 1)/2
 *
 * @param x The element at which the pdf is evaluated
 * @param sigma The standard deviation of the equivalent Gaussian distribution
 * @param q The modulus of the ring Z_q
 * @return The value for the pdf of the uniform distribution at x
 */

double uniform_pdf(long x, double sigma, int n, long q);

/*
 * Translates the vector elements from [a, b] into an index of an array.
 *
 * @param elem The vector to be translated
 * @param q The modulus of the ring Z_q
 * @param a Starting index
 * @param b Ending index
 * @return The index of an array corresponding to the vector elements from [a, b]
 */

size_t index(vec_t elem, long q, int a, int b);

/*
 * Translates an index of an array into vector elements in [a, b].
 *
 * @param res The destination of the vector elements
 * @parma idx The index corresponding to the vector elements to translate
 * @param q The modulus of the ring Z_q
 * @param a Starting index
 * @param b Ending index
 */

void unindex(vec_t res, size_t idx, long q, int a, int b);

/*
 * Checks that the vector has zero elements between a and b.
 *
 * @param elem The vector to be checked
 * @param a Starting index
 * @param b Ending index
 * @return 1 (true) if the vector elements are zero between [a, b], 0 (false)
 *         otherwise
 */

int zero(vec_t elem, int a, int b);

/*
 * Checks that the vector u is the same as the vector v between a and b.
 *
 * @param u The first vector to be checked
 * @param v The second vector to be checked
 * @param a Starting index
 * @param b Ending index
 * @return 1 (true) if the vectors elements are the same between [a, b], 0
 *         (false) otherwise
 */

int equals(vec_t u, vec_t v, int a, int b);

/*
 * Stand-alone C implementation of the error function erf(x)
 *  Source : Abramowitz and Stegun, Handbook of Mathematical Functions, 1965
 *
 * @param x The element at which the erf is evaluated
 * @return The value of the erf at x
 */

double custom_erf(double x);

#endif /* MATHS_H_ */
