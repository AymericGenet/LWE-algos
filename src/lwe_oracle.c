/*
 * lwe_oracle_test.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "lwe_oracle.h"
#include "math.h"
#include <stdio.h>

int lwe_oracle(math_t * res, long * s, int n, long q, double sigma) {
    size_t i;

    res[n].value = 0;
    for (i = 0; i < n; ++i) {

        #ifdef _WIN32
        /* gathers random data */
        if (!CryptGenRandom(hCryptProv, 4, res[i].bytes)) {
             return 0; /* false */
        }
        #endif /* _WIN32 */

        /* TODO: ifdef unix */

        res[i].value %= q;
        res[n].value += (res[i].value * s[i]) % q;
    }
    res[n].value = (res[n].value + NOISE_FUNCTION(sigma)) % q;

    return 1; /* true */
}
