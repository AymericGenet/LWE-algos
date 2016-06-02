/*
 * main.c : typical execution of the BKW algorithm
 *
 *  Created on: May 20, 2016
 *      Author: Aymeric Genet
 */

#include <stdio.h>
#include <stdlib.h>
#include "src/bkw.h"
#include "src/math.h"
#include "src/lwe.h"
#include <math.h>
#include <float.h>
#include <time.h>


#define MAX_RANGE 100

int ns[] = {8, 9, 10, 11, 12, 13, 14, 15, 16};
int qs[] = {67, 83, 101, 127, 149, 197, 227, 257};
int bs[] = {3, 3, 4, 4, 3, 4, 4, 4};
int as[] = {3, 3, 3, 3, 4, 4, 4, 4};

int idx = 0; /* from 0 to 8 */
int m = 32;

int main(int argc, char *argv[]) {
    size_t i, j, k;
    clock_t time;
    int n, b, a, d;
    int successes[] = {0, 0};
    long q;
    unsigned long depth;
    double lf1_sec = 0.0, ll_sec = 0.0, fft_sec = 0.0;
    math_t * guess;
    math_t * vec;
    math_t ** res;
    math_t ** F;
    math_t *** aux;
    math_t *** T;

    init_random();

    /* ================================ INIT ================================ */
    n = ns[idx], q = qs[idx], a = as[idx], b = bs[idx], d = 2;
    /*n = 6, q = 5, a = 2, b = 3, d = 2;*/
    depth = (unsigned long) pow(q, b);
    if (pow(q, b) - depth != 0) {
        fprintf(stderr, "main: depth overflow\n");
        exit(1);
    }

    for (k = 0; k < MAX_RANGE; ++k) {
        T = malloc(a * sizeof(math_t **));
        F = malloc(m * sizeof(math_t *));
        res = malloc(m * sizeof(math_t *));
        guess = calloc((d - 1), sizeof(math_t));
        vec = malloc(n * sizeof(math_t));

        aux = malloc(2 * sizeof(math_t **));
        aux[0] = malloc(a * sizeof(math_t *));
        aux[1] = malloc(2 * sizeof(math_t *));

        for (i = 0; i < a; ++i) {
            aux[0][i] = calloc((n + 1), sizeof(math_t));
            T[i] = calloc(depth, sizeof(math_t *));
            if (T[i] == NULL) {
                fprintf(stderr, "main: could not allocate enough memory\n");
                exit(1);
            }
        }

        for (i = 0; i < m; ++i) {
            F[i] = calloc(d, sizeof(math_t));
            res[i] = calloc((n + 1), sizeof(math_t));
        }

        aux[1][0] = calloc(m, sizeof(math_t));
        aux[1][1] = calloc((d - 1), sizeof(math_t));


        /* ======================== SAMPLE REDUCTION ======================== */

        /* defines secret and sigma in lwe_oracle */
        secret = malloc(n * sizeof(long));
        for (i = 0; i < n; ++i) {
            read_random(vec + i);
            secret[i] = vec[i].value % q;
        }
        sigma = ((double) q)/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

        /* runs algorithm to recover m samples */
        time = clock();
        for (i = 0; i < m; ++i) {
            bkw_lf1(res[i], n, q, b, d - 1, a, T, aux[0]);
        }
        time = clock() - time;
        lf1_sec += time/((double) CLOCKS_PER_SEC);

        /* reduces samples */
        for (i = 0; i < m; ++i) {
            for (j = 0; j < d; ++j) {
                F[i][j].value = res[i][j + n + 1 - d].value;
            }
        }


        /* ======================= HYPOTHESIS TESTING ======================= */

        /* runs hypothesis testing */
        time = clock();
        bkw_hypo_testing(guess, F, d, m, q, sqrt(pow(2, a - 1)) * sigma, aux[1]);
        time = clock() - time;
        ll_sec += time/((double) CLOCKS_PER_SEC);

        /* checks correct answer */
        successes[0]++;
        for (i = 0; i < d - 1; ++i) {
            if (guess[i].value != secret[n - d + 1 + i]) {
                successes[0]--;
                break;
            }
        }

        /* ===================== FAST FOURIER TRANSFORM ===================== */

        /* runs hypothesis testing */
        time = clock();
        bkw_fft(guess, F, d, m, q);
        time = clock() - time;
        fft_sec += time/((double) CLOCKS_PER_SEC);

        /* checks correct answer */
        successes[1]++;
        for (i = 0; i < d - 1; ++i) {
            if (guess[i].value != secret[n - d + 1 + i]) {
                successes[1]--;
                break;
            }
        }

        /* =============================== END ============================== */


        bkw_free_log();

        for (i = 0; i < a; ++i) {
            free(aux[0][i]);
            for (j = 0; j < depth; ++j) {
                free(T[i][j]);
            }
            free(T[i]);
        }

        for (i = 0; i < m; ++i) {
            free(F[i]);
            free(res[i]);
        }

        free(aux[1][0]);
        free(aux[1][1]);
        free(aux[0]);
        free(aux[1]);
        free(aux);

        free(vec);
        free(guess);
        free(res);
        free(F);
        free(T);
        free(secret);
    }

    printf("For m = %i, n = %i, we have :\n\n", m, n);
    printf("\t# of success (Log-likelihood)    : %i / %i\n", successes[0], MAX_RANGE);
    printf("\t# of success (Fourier transform) : %i / %i\n", successes[0], MAX_RANGE);

    printf("\n\taverage time for collecting samples : %f s", lf1_sec / MAX_RANGE);
    printf("\n\taverage time for log-likelihood     : %f s", ll_sec / MAX_RANGE);
    printf("\n\taverage time for Fourier transform  : %f s", fft_sec / MAX_RANGE);

    printf("\n\nPoutcha c'est la meilleure <3\n");
    close_random();
    return 0;
}
