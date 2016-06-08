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


#define MAX_RANGE 2

int ns[] = {6, 7, 8, 9, 10, 11, 12};
int qs[] = {37, 53, 67, 83, 101, 127, 149};
int as[] = {4, 4, 5, 5, 5, 5, 5};

int idx = 0; /* from 0 to 8 */
int m = 5;

int main(int argc, char *argv[]) {
    size_t i, j, k;
    clock_t time;
    int n, b, a, d;
    int successes = 0;
    long q;
    unsigned long depth;
    double lf1_sec = 0.0, sol_sec = 0.0, rdm = 0.0;
    vec_t guess;
    vec_t vec;
    vec_t * res;
    vec_t * F;
    math_t *** aux;
    lwe_t lwe;
    bkw_t bkw;

    init_random();

    /* ================================ INIT ================================ */
    n = ns[idx], q = qs[idx], a = as[idx], b = n/a, d = 1;
    /*n = 6, q = 13, a = 2, b = 3, d = 1;*/
    depth = (unsigned long) pow(q, b);
    if (pow(q, b) - depth != 0) {
        fprintf(stderr, "main: depth overflow\n");
        exit(1);
    }

    for (k = 0; k < MAX_RANGE; ++k) {
        F = malloc(m * sizeof(vec_t));
        res = malloc(m * sizeof(vec_t));
        guess = calloc(d, sizeof(math_t));
        vec = malloc(n * sizeof(math_t));

        aux = malloc(2 * sizeof(math_t **));
        aux[0] = malloc(a * sizeof(math_t *));
        aux[1] = malloc(2 * sizeof(math_t *));

        for (i = 0; i < a; ++i) {
            aux[0][i] = calloc((n + 1), sizeof(math_t));
        }

        for (i = 0; i < m; ++i) {
            F[i] = calloc(d + 1, sizeof(math_t));
            res[i] = calloc((n + 1), sizeof(math_t));
        }

        aux[1][0] = calloc(m, sizeof(math_t));
        aux[1][1] = calloc(d, sizeof(math_t));


        /* ======================== SAMPLE REDUCTION ======================== */

        /* defines secret and sigma in lwe_oracle */
        secret = malloc(n * sizeof(long));
        for (i = 0; i < n; ++i) {
            do {
                read_drandom(&rdm);
                secret[i] = rdm * q;
            } while (secret[i] == 0);
        }
        sigma = ((double) q)/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

        lwe_create(&lwe, n, q, rounded_gaussian, sqrt(pow(2, a)) * (sigma/q));
        bkw_create(&bkw, lwe, a, d, m);

        /* runs algorithm to recover m samples */
        time = clock();
        for (i = 0; i < m; ++i) {
            bkw_lf1(res[i], bkw, a, aux[0]);
        }
        time = clock() - time;

        lf1_sec += time/((double) CLOCKS_PER_SEC);

        /* reduces samples */
        for (i = 0; i < m; ++i) {
            for (j = 0; j < d + 1; ++j) {
                F[i][j] = res[i][j + n - d];
            }
        }


        /* ========================== SOLVING PART ========================== */

        /* runs hypothesis testing or FFT */
        time = clock();
        bkw_hypo_testing(guess, F, bkw, aux[1]);
        time = clock() - time;
        sol_sec += time/((double) CLOCKS_PER_SEC);

        /* checks correct answer */
        successes++;
        for (i = 0; i < d; ++i) {
            if (guess[i] != secret[n - d + i]) {
                successes--;
                break;
            }
        }


        /* =============================== END ============================== */

        bkw_free(&bkw);
        bkw_free_log();

        for (i = 0; i < a; ++i) {
            free(aux[0][i]);
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
        free(secret);
    }

    printf("For m = %i, n = %i, we have :\n\n", m, n);
    printf("\t# of success : %i / %i\n", successes, MAX_RANGE);

    printf("\n\taverage time for collecting samples : %f s", lf1_sec / MAX_RANGE);
    printf("\n\taverage time for solving part       : %f s", sol_sec / MAX_RANGE);

    printf("\n\n");
    close_random();
    return 0;
}
