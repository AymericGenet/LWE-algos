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

int main(int argc, char *argv[]) {
    size_t i, j;
    clock_t time;
    int n, b, a, m, d;
    long q;
    unsigned long depth;
    math_t * guess;
    math_t * vec;
    math_t ** res;
    math_t ** F;
    math_t *** aux;
    math_t *** T;


    /* ================================ INIT ================================ */
    n = 8, q = 67, b = 2, d = 3;
    a = n/b;
    depth = (unsigned long) pow(q, b);
    m = (int) 100;

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
    }

    for (i = 0; i < m; ++i) {
        F[i] = calloc(d, sizeof(math_t));
        res[i] = calloc((n + 1), sizeof(math_t));
    }

    aux[1][0] = calloc(m, sizeof(math_t));
    aux[1][1] = calloc((d - 1), sizeof(math_t));

    init_random();


    /* ========================== SAMPLE REDUCTION ========================== */

    /* defines secret and sigma in lwe_oracle */
    secret = malloc(n * sizeof(long));
    for (i = 0; i < n; ++i) {
        read_random(vec + i);
        secret[i] = vec[i].value % q;
    }
    sigma = ((double) q)/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    /* runs algorithm to recover m samples */
    printf("Collecting %i samples... ", m);
    time = clock();
    for (i = 0; i < m; ++i) {
        bkw_lf1(res[i], n, q, b, a - 1, T, aux[0]);
    }
    time = clock() - time;
    printf("done !\nTime : %f s\n", time/((double) CLOCKS_PER_SEC));

    /* reduces samples */
    for (i = 0; i < m; ++i) {
        for (j = 0; j < d; ++j) {
            F[i][j].value = res[i][j + (a - 1)*b].value;
        }
    }


    /* ========================= HYPOTHESIS TESTING ========================= */

    /* runs hypothesis testing */
    printf("\nRunning log-likelihood hypothesis testing... ");
    time = clock();
    bkw_hypo_testing(guess, F, d, m, q, sqrt(pow(2, a - 1)) * sigma, aux[1]);
    time = clock() - time;
    printf("done !\nTime : %f s\n", time/((double) CLOCKS_PER_SEC));

    /* prints best guess */
    printf("\n\t best guess     v = [ ");
    for (i = 0; i < d - 1; ++i) {
        printf("%lu ", guess[i].value);
    }
    printf("]\n");
    printf("\t correct secret s = [ ");
    for (i = n - d + 1; i < n; ++i) {
        printf("%lu ", secret[i]);
    }
    printf("]\n\n");


    /* ======================= FAST FOURIER TRANSFORM ======================= */

    /* runs hypothesis testing */
    printf("\nRunning FFT hypothesis testing... ");
    time = clock();
    bkw_fft(guess, F, d, m, q);
    time = clock() - time;
    printf("done !\nTime : %f s\n", time/((double) CLOCKS_PER_SEC));

    /* prints best guess */
    printf("\n\t best guess     v = [ ");
    for (i = 0; i < d - 1; ++i) {
        printf("%lu ", guess[i].value);
    }
    printf("]\n");
    printf("\t correct secret s = [ ");
    for (i = n - d + 1; i < n; ++i) {
        printf("%lu ", secret[i]);
    }
    printf("]\n\n");

    /* ================================= END ================================ */

    close_random();
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

    return 0;
}
