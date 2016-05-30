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

int read_sample(math_t * res, FILE * file, long q, int n) {
    size_t i, j;
    char num[128], c;

    memset(num, 0, 128);
    j = 0;
    i = 0;
    do {
        c = fgetc(file);
        switch(c) {
        case ' ':
            num[j] = '\0';
            res[i].value = strtol(num, NULL, 10);
            memset(num, 0, 128);
            i++;
            j = 0;
            break;

        case '\n':
            res[i].value = strtol(num, NULL, 10);
            memset(num, 0, 128);
            return 1; /* true */

        default:
            if (j >= 128) {
                return 0;
            }
            num[j] = c;
            j++;
        }
    } while (c != EOF);

    return 0; /* false */
}

int append_sample(FILE * file, math_t * vec, int n) {
    size_t i;
    int success;

    for (i = 0; i < n; ++i) {
        success = fprintf(file, "%lu%c", vec[i].value, (i == n - 1) ? '\n' : ' ');
    }

    return success;
}

int close_table(FILE ** file) {
    return fclose(*file);
}
