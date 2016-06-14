/*
 * misc.h
 *
 *  Created on: May 29, 2016
 *      Author: Aymeric Genet
 */

#ifndef MISC_H_
#define MISC_H_

#define MAX_SAMPLE_DIGIT 6

#include <stdio.h>
#include "math.h"

/*
 * Initializes the reading of /dev/urandom/.
 *
 * @return 1 (true) if opening /dev/urandom/ was successful, 0 (false) othwerise
 */

int init_random(char * path);

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

int open_table(FILE ** file, char * path);

int read_sample(vec_t res, FILE * file, long q, int n);

int append_sample(FILE * file, vec_t vec, int n);

int close_table(FILE ** file);

#endif /* MISC_H_ */
