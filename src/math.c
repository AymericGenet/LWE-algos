/*
 * math.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "math.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#ifdef _WIN32

HCRYPTPROV hCryptProv = NULL;

int win32_init() {
    if (!CryptAcquireContext(&hCryptProv, NULL, NULL, PROV_RSA_FULL, 0)) {
        return 0; /* false */
    }
    return 1; /* true */
}

int win32_clear() {
    if (!CryptReleaseContext(hCryptProv, 0)) {
        return 0; /* false */
    }
    return 1; /* true */
}

#endif /* _WIN32 */


long discrete_gaussian(double sigma) {
    return 0;
}

long rounded_gaussian(double sigma) {
    return 0;
}

long uniform(double sigma) {
    int beta;
    math_t tmp;

    beta = (sqrt(12*sigma*sigma + 1) - 1)/2;


    #ifdef _WIN32
    /* gathers random data */
    if (!CryptGenRandom(hCryptProv, 4, tmp.bytes)) {
         return 0;
    }
    #endif /* _WIN32 */


    /* TODO: ifdef unix */

    return ((signed long) tmp.value) % beta;
}
