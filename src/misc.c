/*
 * misc.c
 *
 *  Created on: May 29, 2016
 *      Author: Aymeric Genet
 */

#include "misc.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <fcntl.h>
#include <unistd.h>
#include <float.h>
#include <stdio.h>

static int random_file;


int init_random(char * path) {
    random_file = open(path, O_RDONLY);
    return random_file;
}

int read_random(math_t * dest) {
    int success;
    success = read(random_file, dest, sizeof *dest);
    return success;
}

int read_drandom(double * dest) {
    unsigned long random;
    int success;

    success = read(random_file, &random, sizeof(unsigned long));
    if (success < 0) {
        return success;
    }

    *dest = random/((double) ULONG_MAX);

    return success;
}

int close_random() {
    return close(random_file);
}

int open_table(FILE ** file, char * path) {
    (*file) = fopen(path, "a+");
    return (*file) == NULL;
}

int read_sample(vec_t res, FILE * file, long q, int n) {
    size_t i, j;
    char num[MAX_SAMPLE_DIGIT + 1], c;

    /* initializes memory to 0 */
    memset(num, 0, MAX_SAMPLE_DIGIT + 1);
    j = 0;
    i = 0;

    /* reads one line character by character */
    do {
        c = fgetc(file);
        switch(c) {
        case ' ':
            /* halts parsing and converts into long */
            num[j] = '\0';
            res[i] = strtol(num, NULL, 10);

            /* resets memory to 0 */
            memset(num, 0, MAX_SAMPLE_DIGIT + 1);
            j = 0;

            /* continues reading */
            i++;
            break;

        case '\n':
            /* converts to long */
            res[i] = strtol(num, NULL, 10);

            /* resets memory to 0 */
            memset(num, 0, MAX_SAMPLE_DIGIT + 1);
            return 1; /* true */

        default:
            /* aborts if too many digits */
            if (j > MAX_SAMPLE_DIGIT) {
                return 0;
            }

            /* appends number */
            num[j] = c;
            j++;
        }
    } while (c != EOF);

    return 0; /* false */
}

int append_sample(FILE * file, vec_t vec, int n) {
    size_t i;
    int success;

    success = 1; /* true */
    for (i = 0; success && i < n; ++i) {
        success = fprintf(file, "%lu%c", vec[i], (i == n - 1) ? '\n' : ' ');
    }

    return success;
}

int close_table(FILE ** file) {
    return fclose(*file);
}
