/*
 * lwe_oracle_test.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "lwe_oracle.h"
#include "math.h"

long * secret;
double sigma;

int lwe_oracle_predef(math_t * res, long * s, int n, long q, double sig) {
    size_t i;

    res[n].value = NOISE_FUNCTION(sig);
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

    return 1; /* true */
}

int lwe_oracle(math_t * res, int n, long q) {
    return lwe_oracle_predef(res, secret, n, q, sigma);
}
