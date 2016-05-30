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

/* Represents an integer either as an array of bytes or a long value */
/* FIXME : obsolete, change into unsigned long */
typedef union {
    unsigned char bytes[sizeof(long)];
    unsigned long value;
} math_t;

/*
 * Initializes the reading of /dev/urandom/.
 *
 * @return 1 (true) if opening /dev/urandom/ was successful, 0 (false) othwerise
 */

int init_random();

/*
 * Reads 4 random bytes from /dev/urandom/ and writes it in dest.
 *  WARNING : must call init_random() beforehand, or the function blocks.
 *
 * @param dest The destination for the random bytes
 * @return 1 (true) if reading /dev/urandom/ was successful, 0 (false) otherwise
 */

int read_random(math_t * dest);

/*
 * Reads 4 random bytes from /dev/urandom/ and writes a random floating between
 * [0, 1] in dest.
 *  WARNING : must call init_random() beforehand, or the function blocks.
 *
 * @param dest The destination for the random double
 * @return 1 (true) if reading /dev/urandom/ was successful, 0 (false) otherwise
 */

int read_drandom(double * dest);

/*
 * Closes the reading of /dev/urandom/.
 *
 * @return 1 (true) if closing /dev/urandom/ was successful, 0 (false) othwerise
 */

int close_random();

/*
 * Samples a number from the discrete Gaussian distribution.
 *
 * @param sigma The standard deviation of the Gaussian distribution
 * @param q The modulus of the ring Z_q
 * @output A sample from the discrete Gaussian distribution
 */

long discrete_gaussian(double sigma, long q);

/*
 * Computes the probability density function of the discrete Gaussian
 * distribution.
 *
 * @param x The element at which the pdf is evaluated
 * @param sigma The standard deviation of the Gaussian distribution
 * @param q The modulus of the ring Z_q
 * @return The value for the cdf of the rounded Gaussian at x
 */

double discrete_gaussian_pdf(long x, double sigma, long q);

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

double rounded_gaussian_pdf(long x, double sigma, long q, long k);

/*
 * Samples a number from the rounded Gaussian distribution.
 *
 * @param sigma The standard deviation of the Gaussian distribution
 * @param q The modulus of the ring Z_q
 * @output A sample from the discrete Gaussian distribution
 */

long rounded_gaussian(double sigma, long q);

/*
 * Samples a number from the uniform distribution between [-beta, beta], where
 *      beta = (sqrt(12*sigma*sigma + 1) - 1)/2
 *
 * @param sigma The standard deviation of the equivalent Gaussian distribution
 * @param q The modulus of the ring Z_q
 * @output A sample from the discrete uniform distribution
 */

long uniform(double sigma, long q);

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

double uniform_pdf(long x, double sigma, long q);

/*
 * Translates the vector elements from [a, b] into an index of an array.
 *
 * @param elem The vector to be translated
 * @param q The modulus of the ring Z_q
 * @param a Starting index
 * @param b Ending index
 * @return The index of an array corresponding to the vector elements from [a, b]
 */

size_t index(math_t * elem, long q, int a, int b);

/*
 * Translates an index of an array into vector elements in [a, b].
 *
 * @param res The destination of the vector elements
 * @parma idx The index corresponding to the vector elements to translate
 * @param q The modulus of the ring Z_q
 * @param a Starting index
 * @param b Ending index
 */

void unindex(math_t * res, size_t idx, long q, int a, int b);

/*
 * Checks that the vector has zero elements between a and b.
 *
 * @param elem The vector to be checked
 * @param a Starting index
 * @param b Ending index
 * @return 1 (true) if the vector elements are zero between [a, b], 0 (false)
 *         otherwise
 */

int zero(math_t * elem, int a, int b);

int equals(math_t * u, math_t * v, int a, int b);

/*
 * Stand-alone C implementation of the error function erf(x)
 *  Source : Abramowitz and Stegun, Handbook of Mathematical Functions, 1965
 *
 * @param x The element at which the erf is evaluated
 * @return The value of the erf at x
 */

double custom_erf(double x);

#endif /* MATHS_H_ */
