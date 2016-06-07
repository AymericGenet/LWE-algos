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


int open_table(FILE ** file, char * path) {
    (*file) = fopen(path, "a+");
    return (*file) == NULL;
}

int read_sample(vec_t res, FILE * file, long q, int n) {
    size_t i, j;
    char num[MAX_SAMPLE_DIGIT + 1], c;

    memset(num, 0, MAX_SAMPLE_DIGIT + 1);
    j = 0;
    i = 0;
    do {
        c = fgetc(file);
        switch(c) {
        case ' ':
            num[j] = '\0';
            res[i] = strtol(num, NULL, 10);
            memset(num, 0, MAX_SAMPLE_DIGIT + 1);
            i++;
            j = 0;
            break;

        case '\n':
            res[i] = strtol(num, NULL, 10);
            memset(num, 0, MAX_SAMPLE_DIGIT + 1);
            return 1; /* true */

        default:
            if (j > MAX_SAMPLE_DIGIT) {
                return 0;
            }
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
