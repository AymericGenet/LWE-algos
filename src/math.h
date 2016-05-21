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
 */

int init_random();

/*
 * Reads 4 random bytes from /dev/urandom/ and writes it in dest.
 *  WARNING : must call init_random() beforehand, or the function blocks.
 */

int read_random(math_t * dest);

/*
 * Reads 4 random bytes from /dev/urandom/ and writes a random floating between
 * [0, 1] in dest.
 *  WARNING : must call init_random() beforehand, or the function blocks.
 */

int read_drandom(double * dest);

/*
 * Closes the reading of /dev/urandom/.
 */

int close_random();

/*
 * Samples a number from the discrete Gaussian distribution.
 */

long discrete_gaussian(double sigma);

/*
 * Computes the cumulative distribution function of the rounded Gaussian
 * distribution.
 */

double rounded_gaussian_cdf(double x, double sigma);

/*
 * Computes the probability density function of the rounded Gaussian
 * distribution.
 */

double rounded_gaussian_pdf(long x, double sigma, long q, long k);

/*
 * Samples a number from the rounded Gaussian distribution.
 */

long rounded_gaussian(double sigma, long q);

/*
 * Samples a number from the uniform distribution between [-beta, beta], where
 *      beta = (sqrt(12*sigma*sigma + 1) - 1)/2
 */

long uniform(double sigma, long q);

/*
 * Computes the probability density function of the uniform distribution between
 * [-beta, beta], where
 *      beta = (sqrt(12*sigma*sigma + 1) - 1)/2
 */

double uniform_pdf(long x, double sigma, long q);

/*
 * Translates a vector into an index of an array.
 */

size_t index(math_t * elem, long q, int a, int b);

/*
 * Translates an index of an array into a vector.
 */

void unindex(math_t * res, size_t idx, long q, int a, int b);

/*
 * Returns 1 (true) if the vector has zero elements between a and b, 0 (false)
 * otherwise.
 */

int zero(math_t * elem, int a, int b);

/*
 * Stand-alone C implementation of the error function erf(x)
 *  Source : Abramowitz and Stegun, Handbook of Mathematical Functions, 1965
 */

double custom_erf(double x);

#endif /* MATHS_H_ */
