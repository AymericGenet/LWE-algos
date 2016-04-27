/*
 * maths.h : Definitions for mathematical operations.
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#ifndef MATH_H_
#define MATH_H_


#ifdef _WIN32

#include <windows.h>
#include <wincrypt.h>

#endif /* _WIN32 */


#ifdef __unix__

/* TODO: define /dev/urandom */

#endif /* __unix__ */


#define PI_VAL 3.14159265358979323846

/* Represents an integer either as an array of bytes or a long value */
typedef union {
    unsigned char bytes[4];
    unsigned long value;
} math_t;


#ifdef _WIN32

/* Cryptographic context provider for using wincrypt.h */
extern HCRYPTPROV hCryptProv;

/*
 * Initializes cryptographic context provider, to be called before any operation
 * when using Wincrypt API.
 */

int win32_init();

/*
 * Clears cryptographic context provider, to be called when finished using the
 * Wincrypt API.
 */

int win32_clear();

#endif /* _WIN32 */


/*
 * To be defined...
 */

long discrete_gaussian(double sigma);

/*
 * To be defined...
 */

long rounded_gaussian(double sigma);

/*
 * Generates number uniformly in [-beta, beta], where
 *      beta = (sqrt(12*sigma*sigma + 1) - 1)/2
 */

long uniform(double sigma);

#endif /* MATHS_H_ */
