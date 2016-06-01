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

int open_table(FILE ** file, char * path);

int read_sample(math_t * res, FILE * file, long q, int n);

int append_sample(FILE * file, math_t * vec, int n);

int close_table(FILE ** file);

#endif /* MISC_H_ */
