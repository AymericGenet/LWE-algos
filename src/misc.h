/*
 * misc.h
 *
 *  Created on: May 29, 2016
 *      Author: Aymeric Genet
 */

#ifndef MISC_H_
#define MISC_H_

#define MAX_SAMPLE_DIGIT 6

#include "math.h"

#include <stdio.h>

/*
 * Initializes the reading of a file containing random elements. Typical uses
 * include "/dev/random" and "/dev/urandom".
 *
 * @return 1 (true) if opening the file was successful, 0 (false) othwerise
 */

int init_random(char * path);

/*
 * Reads 4 random bytes from the file initiated by init_random() and writes it
 * in dest.
 *  WARNING : must call init_random() beforehand, or the function blocks.
 *
 * @param dest The destination for the random bytes
 * @return 1 (true) if reading the file was successful, 0 (false) otherwise
 */

int read_random(math_t * dest);

/*
 * Reads 4 random bytes from the file initiated by init_random() and writes a
 * random floating between [0, 1] in dest.
 *  WARNING : must call init_random() beforehand, or the function blocks.
 *
 * @param dest The destination for the random double
 * @return 1 (true) if reading the file was successful, 0 (false) otherwise
 */

int read_drandom(double * dest);

/*
 * Closes the reading of the file initiated by init_random().
 *
 * @return 1 (true) if closing the file was successful, 0 (false) othwerise
 */

int close_random();

/*
 * Opens a persistent table T^l containing samples for sample reduction.
 *
 * @param file File to be read after opening
 * @param path Path to a persistent table T^l
 * @return 1 (true) if opening the file was successful, 0 (false) otherwise
 */

int open_table(FILE ** file, char * path);

/*
 * Reads next sample from file and puts it into res.
 *
 * @param res Vector containing the sample read
 * @param file File to be read
 * @return 1 (true) if reading was successful, 0 (false) otherwise
 */

int read_sample(vec_t res, FILE * file, long q, int n);

/*
 * Appends a sample in a file.
 *
 * @param file File in which sample has to be appended
 * @param vec Sample to be appended
 * @param n Size of the sample
 * @return 1 (true) if appending the sample was successful, 0 (false) otherwise
 */

int append_sample(FILE * file, vec_t vec, int n);

/*
 * Closes the file to the persistent table T^l.
 *
 * @param file File to close
 * @return 1 (true) if closing the file was successful, 0 (false) othwerise
 */

int close_table(FILE ** file);

#endif /* MISC_H_ */
