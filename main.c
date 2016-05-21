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
#include "src/lwe_oracle.h"
#include <math.h>
#include <float.h>

int main(int argc, char *argv[]) {
    size_t i, j;
    int n, b, a, m, d, k;
    long q, depth, basis, noise;
    math_t * guess;
    math_t ** res;
    math_t ** F;
    math_t *** aux;
    math_t *** T;
    double ** S;
    double best, sum;
    unsigned int * stats;


    /* ================================ INIT ================================ */
    n = 6, q = 13, b = 3, d = 4;
    a = n/b;
    depth = 2197;
    m = 10;

    T = malloc(a * sizeof(math_t **));
    F = malloc(m * sizeof(math_t *));
    S = malloc(m * sizeof(double *));
    res = malloc(m * sizeof(math_t *));
    stats = malloc(q * sizeof(unsigned long));
    guess = malloc((d - 1) * sizeof(math_t));

    aux = malloc(2 * sizeof(math_t **));
    aux[0] = malloc(a * sizeof(math_t *));
    aux[1] = malloc(2 * sizeof(math_t *));

    for (i = 0; i < a; ++i) {
        aux[0][i] = malloc((n + 1) * sizeof(math_t));
        T[i] = malloc(depth * sizeof(math_t *));

        for (j = 0; j < depth; ++j) {
            T[i][j] = NULL;
        }
        for (j = 0; j < n + 1; ++j) {
            aux[0][i][j].value = 0;
        }
    }

    for (i = 0; i < m; ++i) {
        F[i] = malloc(d * sizeof(math_t));
        S[i] = malloc(depth * sizeof(double));
        res[i] = malloc((n + 1) * sizeof(math_t));

        for (j = 0; j < depth; ++j) {
            S[i][j] = 0.0;
        }
    }

    aux[1][0] = malloc(m * sizeof(math_t));
    aux[1][1] = malloc((d - 1) * sizeof(math_t));
    for (i = 0; i < d - 1; ++i) {
        aux[1][1][i].value = 0;
        guess[i].value = 0;
    }

    init_random();


    /* ========================== SAMPLE REDUCTION ========================== */

    /* defines secret and sigma in lwe_oracle */
    secret = malloc(n * sizeof(long));
    secret[0] = 1;
    secret[1] = 2;
    secret[2] = 3;
    secret[3] = 4;
    secret[4] = 5;
    secret[5] = 6;
    sigma = ((double) q)/(sqrt(2 * PI_VAL * n) * log(n) * log(n));

    /* runs algorithm to recover m samples */
    printf("Collecting %i samples... ", m);
    for (i = 0; i < m; ++i) {
        bkw_lf1(res[i], n, q, b, a - 1, T, aux[0]);
    }
    printf("done !\n");

    /* reduces samples */
    for (i = 0; i < m; ++i) {
        for (j = 0; j < d; ++j) {
            F[i][j].value = res[i][j + (a - 1)*b].value;
        }
    }

    /* prints distribution */
    printf("\nAimed distribution :\n");
    for (i = 0; i < q; ++i) {
        stats[i] = 0;
    }
    for (i = 0; i < 1000; ++i) {
        stats[rounded_gaussian(sqrt(pow(2, a - 1)) * sigma, q)]++;
    }
    for (i = 0; i < q; ++i) {
        sum = (stats[i] / 10.0);
        for (j = 0; j < sum; ++j) {
            printf("=");
        }
        printf(" (%f)\n", rounded_gaussian_pdf(i, sqrt(pow(2, a - 1)) * sigma, q, 10));
    }

    /* prints noise extraction distribution */
    printf("\nObtained distribution :\n");
    for (i = 0; i < q; ++i) {
        stats[i] = 0;
    }
    for (i = 0; i < m; ++i) {
        noise = 0;
        for (j = 0; j < d - 1; ++j) {
            noise = (noise + secret[n - d + j + 1]*F[i][j].value) % q;
        }
        noise = (F[i][d - 1].value + q - noise) % q;
        stats[noise]++;
    }
    for (i = 0; i < q; ++i) {
        sum = 100.0 * (stats[i] / (double) m);
        for (j = 0; j < sum; ++j) {
            printf("=");
        }
        printf(" (%f)\n", stats[i] / (double) m);
    }


    /* ========================= HYPOTHESIS TESTING ========================= */

    /* runs hypothesis testing */
    printf("\nRunning hypothesis testing... ");
    bkw_hypo_testing(S, F, d, m, q, sqrt(pow(2, a - 1)) * sigma, aux[1]);
    printf("done\n");

    /* recovers the best vector according to its score */
    best = -DBL_MAX;
    k = -1;
    for (i = 1; i < depth; ++i) {
        sum = 0.0;
        for (j = 0; j < m; ++j) {
            sum += S[j][i];
        }
        if (sum >= best) {
            k = i;
            best = sum;
        }
    }
    unindex(guess, k, q, n - d + 1, n);

    printf("\n\t best guess     v = [ ");
    for (i = 0; i < d - 1; ++i) {
        printf("%lu ", guess[i].value);
    }
    printf("] (idx %i) score of %.5f\n", k, best);

    /* compares with the correct vector */
    k = 0;
    basis = 1;
    for (j = n - b; j < n; ++j) {
        k += secret[j]*basis;
        basis *= q;
    }
    sum = 0;
    for (j = 0; j < m; ++j) {
        sum += S[j][k];
    }
    printf("\t correct secret s = [ ");
    for (i = n - d + 1; i < n; ++i) {
        printf("%lu ", secret[i]);
    }
    printf("] (idx %i) score of %.5f\n\n", k, sum);


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
        free(S[i]);
        free(res[i]);
    }

    free(aux[1][0]);
    free(aux[1][1]);
    free(aux[0]);
    free(aux[1]);
    free(aux);

    free(guess);
    free(stats);
    free(S);
    free(res);
    free(F);
    free(T);
    free(secret);

    return 0;
}
