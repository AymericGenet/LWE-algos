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
#include "src/misc.h"
#include "src/lwe.h"
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>


int ns[] = {6, 7, 8, 9, 10, 11, 12};
int qs[] = {37, 53, 67, 83, 101, 127, 149};
int as[] = {5, 6, 7, 8, 9, 9, 10};
int bs[] = {1, 1, 1, 1, 1, 1, 1};

int idx = 0; /* from 0 to 8 */
int m = 15;
int MAX_RANGE = 1000;

int main(int argc, char *argv[]) {
    int * stats;
    long noise;
    int i, j, k;
    clock_t time;
    int n, b, a, d;
    int successes = 0;
    long q;
    unsigned long depth = 0, mem_count = 0;
    double avg_oracle_calls = 0.0, avg_memory = 0.0;
    double lf1_sec = 0.0, sol_sec = 0.0, rdm = 0.0;
    double sigma;
    distribution_t distrib;
    vec_t guess;
    vec_t sec;
    vec_t * res;
    vec_t * F;
    math_t *** aux;
    lwe_t * lwe;
    bkw_t * bkw;
    node_t * curr;
    char c;
    int verbose = 0;
    int failure = 0;
    int sample_reduc = 0, hypo_test = 0;

    /* ================================ DATA ================================ */
    distrib = rounded_gaussian;
    while ((c = getopt(argc, argv, "i:n:q:a:m:r:d:h:s:v")) != -1)
        switch (c) {
        case 'v':
            verbose = 1;
            break;
        case 's':
            sample_reduc = atoi(optarg);
            break;
        case 'h':
            hypo_test = atoi(optarg);
            break;
        case 'i':
            idx = atoi(optarg);
            n = ns[idx], q = qs[idx], a = as[idx], b = bs[idx], d = n - (a - 1)*b;
            break;
        case 'n':
            n = atoi(optarg);
            break;
        case 'q':
            q = atoi(optarg);
            break;
        case 'a':
            a = atoi(optarg);
            break;
        case 'm':
            m = atoi(optarg);
            break;
        case 'r':
            MAX_RANGE = atoi(optarg);
            break;
        case 'd':
            switch (atoi(optarg)) {
            case 2:
                distrib = uniform;
                break;

            case 1:
                distrib = discrete_gaussian;
                break;

            case 0:
            default:
                distrib = rounded_gaussian;
            }
            break;
        default:
            exit(1);
    }

    /* checks feasibility of the instance */
    depth = (unsigned long) pow(q, b);
    if (pow(q, b) - depth != 0) {
        fprintf(stderr, "main: depth overflow\n");
        exit(1);
    }
    init_random("/dev/urandom");

    for (k = 0; k < MAX_RANGE; ++k) {
        /* ============================== INIT ============================== */
        F = malloc(m * sizeof(vec_t));
        res = malloc(m * sizeof(vec_t));
        guess = calloc(d, sizeof(math_t));
        
        stats = calloc(q, sizeof(int));
        aux = malloc(2 * sizeof(math_t **));
        aux[0] = malloc(a * sizeof(math_t *));
        aux[1] = malloc(2 * sizeof(math_t *));

        /* defines secret and sigma in lwe_oracle */
        sec = malloc(n * sizeof(math_t));
        for (i = 0; i < n; ++i) {
            do {
                read_drandom(&rdm);
                sec[i] = rdm * (q - 1);
            } while (sec[i] == 0);
        }
        sigma = ((double) q)/(sqrt(2 * PI_VAL * n) * log(n)/log(2) * log(n)/log(2));

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

        /* creates lwe and bkw struct */
        lwe = malloc(sizeof(lwe_t));
        lwe_create(lwe, n, q, distrib, sigma, sec);
        bkw = malloc(sizeof(bkw_t));
        bkw_create(bkw, lwe, a, b, d, m);

        /* runs BKW algorithm to recover m samples */
        lwe_oracle_calls = 0;
        time = clock();
        for (i = 0; i < m; ++i) {
            do {
                switch (sample_reduc) {
                    case 1:
                    bkw_lf2(res[i], bkw, a - 1, aux[0]);
                    break;

                    case 0:
                    default:
                    bkw_lf1(res[i], bkw, a - 1, aux[0]);
                }
            } while (zero(res[i], 0, n));
        }
        time = clock() - time;

        /* computes time and oracle calls average */
        avg_oracle_calls = avg_oracle_calls * k + lwe_oracle_calls;
        avg_oracle_calls = avg_oracle_calls / (k + 1);
        lf1_sec += time/((double) CLOCKS_PER_SEC);

        /* computes required amount of memory */
        mem_count = 0;
        curr = bkw->tab->first;

        do {
            for (i = 0; i < a; ++i) {
                for (j = 0; j < depth; ++j) {
                    if (curr->T[i][j] != NULL) {
                        mem_count++;
                    }
                }
            }
            curr = curr->next;
        } while(curr != NULL);
        avg_memory = avg_memory * k + mem_count;
        avg_memory = avg_memory / (k + 1);

        /* reduces samples */
        for (i = 0; i < m; ++i) {
            for (j = 0; j < d + 1; ++j) {
                F[i][j] = res[i][j + n - d];
            }
        }


        /* ========================== SOLVING PART ========================== */
        /* runs hypothesis testing or FFT */
        time = clock();
        switch (hypo_test) {
            case 1:
            bkw_fft(guess, F, bkw);
            break;

            case 0:
            default:
            bkw_hypo_testing(guess, F, bkw, aux[1]);
        }
        time = clock() - time;
        sol_sec += time/((double) CLOCKS_PER_SEC);

        /* checks correct answer */
        successes++;
        failure = 0;
        for (i = 0; i < d; ++i) {
            if (guess[i] % q != sec[n - d + i] % q) {
                successes--;
                failure = 1;
                break;
            }
        }

        /* on failure, when verbose, prints distributions */
        if (verbose && failure) {
            for (i = 0; i < m; ++i) {
                noise = res[i][n];
                for (j = 0; j < n; ++j) {
                    noise = (q*q + noise - res[i][j]*sec[j]) % q;
                }
                stats[noise]++;
            }

            /* prints secret */
            printf("s = [ ");
            for (i = 0; i < d; ++i) {
                printf("%lu ", sec[n - d + i]);
            }
            printf("] noise : \n");

            for (i = 0; i < q; ++i) {
                printf("%i : \t", i);
                for (j = 0; j < 100.0*stats[i]/((double)m); ++j) {
                    printf("=");
                }
                printf(" (%f vs %f)\n", stats[i]/((double) m), distributions[distrib](i, sigma, a, q));
                stats[i] = 0;
            }
            printf("\n");

            for (i = 0; i < m; ++i) {
                noise = res[i][n];
                for (j = n - d; j < n; ++j) {
                    noise = (q*q + noise - res[i][j]*guess[j - (n - d)]) % q;
                }
                stats[noise]++;
            }

            /* prints guess */
            printf("v = [ ");
            for (i = 0; i < d; ++i) {
                printf("%lu ", guess[i]);
            }
            printf("] noise : \n");

            for (i = 0; i < q; ++i) {
                printf("%i : \t", i);
                for (j = 0; j < 100.0*stats[i]/((double)m); ++j) {
                    printf("=");
                }
                printf(" (%f vs %f)\n", stats[i]/((double) m), distributions[distrib](i, sigma, a, q));
                stats[i] = 0;
            }
            printf("\n");
        }


        /* =============================== END ============================== */
        bkw_free(bkw);
        lwe_free(lwe);
        free(bkw);
        free(lwe);
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

        free(stats);
        free(guess);
        free(res);
        free(sec);
        free(F);
    }

    printf("For m = %i, n = %i, a = %i, we have :\n\n", m, n, a);
    printf("\t# of successes : %i / %i\n", successes, MAX_RANGE);

    printf("\n\taverage time for %s   : %f s", (sample_reduc == 1) ? "LF2" : "LF1", lf1_sec / MAX_RANGE);
    printf("\n\taverage time for %s   : %f s", (hypo_test == 1) ? "FFT" : "L-L", sol_sec / MAX_RANGE);

    printf("\n");

    printf("\n\taverage number of LWE oracle calls    : %f calls", avg_oracle_calls);
    printf("\n\taverage memory for collecting samples : %f samples", avg_memory);

    printf("\n\n");
    close_random();
    return 0;
}
