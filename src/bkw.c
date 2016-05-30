/*
 * bkw.c
 *
 *  Created on: Apr 28, 2016
 *      Author: Aymeric Genet
 */

#include "bkw.h"
#include "math.h"
#include "lwe.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


static double * log_j = NULL;

void bkw_create_table(table_t * tab, int n, long q, int b, int d) {
    size_t i;

    tab->n = n;
    tab->q = q;
    tab->b = b;
    tab->d = d;
    tab->a = n/b;

    tab->first = malloc(sizeof(node_t));
    tab->states = malloc(tab->a * sizeof(int));
    tab->sample = malloc(tab->a * sizeof(math_t *));

    for (i = 0; i < tab->a; ++i) {
        tab->states[i] = 0;
        tab->sample[i] = NULL;
    }

    bkw_create_node(tab->first, tab->a, tab->q, tab->b);
}

void bkw_free_table(table_t * tab) {
    size_t i;

    for (i = 0; i < tab->a; ++i) {
        free(tab->sample[i]);
    }

    free(tab->sample);
    free(tab->states);
    bkw_free_node(tab->first, tab->a, tab->q, tab->b);
    free(tab);
}

void bkw_create_node(node_t * node, int a, long q, int b) {
    size_t i, j;
    long depth;

    depth = pow(q, b);

    node->next = NULL;
    node->T = malloc(a * sizeof(math_t **));
    for (i = 0; i < a; ++i) {
        node->T[i] = malloc(depth * sizeof(math_t *));
        for (j = 0; j < depth; ++j) {
            node->T[i][j] = NULL;
        }
    }
}

void bkw_free_node(node_t * node, int a, int q, int b) {
    size_t i, j;
    long depth;

    depth = pow(q, b);

    if (node->next != NULL) {
        bkw_free_node(node->next, a, q, b);
    }

    for (i = 0; i < a; ++i) {
        for (j = 0; j < depth; ++j) {
            free(node->T[i][j]);
        }
        free(node->T[i]);
    }
    free(node->T);

    free(node);
}

int bkw_algo(math_t * res, int n, int b, int l) {
    return 1;
}

int bkw_lf1(math_t * res, int n, long q, int b, int l, math_t *** T,
            math_t ** aux) {
    size_t i, idx;

    while (1) {
        /* samples from previous layer */
        if (l == 0) {
            return lwe_oracle(res, n, q);
        } else {
            if (!bkw_lf1(aux[0], n, q, b, l - 1, T, aux + 1)) {
                return 0; /* false */
            }
        }

        /* if first elements already 0, returns it */
        if (zero(aux[0], (l - 1) * b, l * b)) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = aux[0][i].value;
            }
            return 1; /* true */
        }

        /* checks if collision */
        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (aux[0][i].value + q - T[l][idx][i].value) % q;
            }
            return 1; /* true */
        }

        /* negates it */
        for (i = 0; i < n + 1; ++i) {
            aux[0][i].value = (q - aux[0][i].value) % q;
        }

        /* checks if collision */
        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (aux[0][i].value + q - T[l][idx][i].value) % q;
            }
            return 1; /* true */
        }

        /* otherwise, stores it, and repeat */
        T[l][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l][idx][i].value = aux[0][i].value;
        }
    }
    return 0; /* should not happen */
}

int bkw_lf2(math_t * res, int n, long q, int b, int l, table_t * tab,
            math_t ** aux) {
    size_t i, idx, n_idx;
    int layer;
    node_t * current;
    node_t * next;
    math_t *** T;
    math_t * sample;

    /* if call to lwe_oracle */
    if (l == 0) {
        return lwe_oracle(res, n, q);
    }

    /* finds layer at which the last sample was output */
    sample = tab->sample[l];
    current = tab->first;
    next = current;
    if (tab->states[l] != -1) {
        layer = 0;
        while (layer != tab->states[l]) {
            current = current->next;
            layer++;
        }
        next = current->next;
    }

    /* if next layer exists, checks if sample still collides there */
    if (sample != NULL && next != NULL) {
        T = next->T;
        (tab->states[l])++;

        /* checks if collision */
        idx = index(sample, q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (sample[i].value + q - T[l][idx][i].value) % q;
            }
            return 1; /* true */
        }

        /* negates it */
        for (i = 0; i < n + 1; ++i) {
            aux[0][i].value = (q - sample[i].value) % q;
        }

        /* checks if collision */
        n_idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][n_idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (aux[0][i].value + q - T[l][n_idx][i].value) % q;
            }
            return 1; /* true */
        }

        /* otherwise, stores it */
        T[l][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l][idx][i].value = sample[i].value;
        }
    }
    /* otherwise, if not first call, creates new layer and puts sample in it */
    else if (sample != NULL) {
        current->next = malloc(sizeof(node_t));
        bkw_create_node(current->next, n/b, q, b);
        T = current->next->T;

        /* stores sample */
        idx = index(sample, q, (l - 1) * b, l * b);
        T[l][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l][idx][i].value = sample[i].value;
        }
    }
    /* upon first call */
    else {
        tab->sample[l] = malloc((n + 1) * sizeof(math_t));
        sample = tab->sample[l];
    }

    T = tab->first->T;
    tab->states[l] = 0;
    while (1) {
        /* samples from previous layer */
        if (!bkw_lf2(aux[0], n, q, b, l - 1, tab, aux + 1)) {
            return 0; /* false */
        }

        /* if first elements already 0, returns it */
        if (zero(aux[0], (l - 1) * b, l * b)) {
            printf("\tl = %i [ ", l);
            for (i = 0; i < n + 1; ++i) {
                printf("%lu ", aux[0][i].value);
                res[i].value = aux[0][i].value;
                sample[i].value = aux[0][i].value;
            }
            printf("]\n");
            tab->states[l] = -1;
            return 1; /* true */
        }

        /* checks if collision */
        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (aux[0][i].value + q - T[l][idx][i].value) % q;
                sample[i].value = aux[0][i].value;
            }
            return 1; /* true */
        }

        /* negates it */
        for (i = 0; i < n + 1; ++i) {
            aux[0][i].value = (q - aux[0][i].value) % q;
        }

        /* checks if collision */
        n_idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][n_idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i].value = (aux[0][i].value + q - T[l][n_idx][i].value) % q;
                sample[i].value = aux[0][i].value;
            }
            return 1; /* true */
        }

        /* otherwise, stores it, and repeat */
        T[l][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l][idx][i].value = (q - aux[0][i].value) % q;
        }
    }
    return 0;
}

void bkw_hypo_testing(math_t * v, math_t ** F, int d, int m, long q,
                      double sigma, math_t ** aux) {
    int i;
    size_t j;
    double p, best, current;

    /* precomputes each value of log(j) */
    if (log_j == NULL) {
        log_j = malloc(q * sizeof(double));
        for (i = 0; i < q; ++i) {
            p = rounded_gaussian_pdf(i, sigma, q, 10*d);
            log_j[i] = log(p) - log((pow(q, d-1) - p)/(double)(pow(q, d) - 1.0));
        }
    }

    /* aux[0][i] = first coordinate of a_i */
    for (i = 0; i < m; ++i) {
        aux[0][i].value = F[i][0].value;
    }

    best = 0.0;
    /* foreach v (aux[1]) in Z_q^d, starting at v = (1 0 0 0...) */
    for (aux[1][0].value = 1; aux[1][d - 2].value != q; ) {
        /* foreach sample, update S[v] with a_i - c_i */
        current = 0.0;
        for (i = 0; i < m; ++i) {
            /* S[v]_i = W(aux[0][i] - c_i) */
            current += log_j[(aux[0][i].value - F[i][d - 1].value + q) % q];
            /* aux[0][i] += first coordinate of a_i */
            aux[0][i].value = (aux[0][i].value + F[i][0].value) % q;
        }
        if (current > best) {
            for (i = 0; i < d - 1; ++i) {
                v[i].value = aux[1][i].value;
            }
            best = current;
        }

        /* increases v by (1 0 0 0...) */
        aux[1][0].value += 1;

        /* when v becomes (q 0 0 0...), increases the next coordinates,
        repeats the process if next coordinates also hit q (max). */
        for (i = 0; aux[1][i].value >= q && i < d - 2 ; ++i) {
            aux[1][i].value = 0;
            /* if there is one coordinate after this one */
            if (i + 1 < d - 1) {
                aux[1][i + 1].value += 1;
                for (j = 0; j < m; ++j) {
                    aux[0][j].value = aux[0][j].value + F[j][i + 1].value;
                }
            } else {
                aux[1][i].value = q; /* to get outside of the loop */
            }
        }
    }
}

void bkw_fft(math_t * v, math_t ** F, int d, int m, long q) {
    size_t i, j, idx, size;
    double best;
    int * n;
    fftw_complex * in;
    fftw_complex * out;
    fftw_plan plan_forward;

    /* creates the input array */
    n = malloc((d - 1) * sizeof(int));
    size = 1;

    for (i = 0; i < d - 1; ++i) {
        size *= q;
        n[i] = q;
    }
    in = fftw_malloc(sizeof(fftw_complex) * size );

    for (i = 0; i < size; ++i) {
        in[i][0] = 0.0;
        in[i][1] = 0.0;
    }

    /* puts value to be transformed */
    for (i = 0; i < m; ++i) {
        idx = 0;
        for (j = 0; j < d - 1; ++j) {
            idx = F[i][j].value + q*idx;
        }
        in[idx][0] += cos(2 * PI_VAL * F[i][d - 1].value /((double) q));
        in[idx][1] += sin(2 * PI_VAL * F[i][d - 1].value /((double) q));
    }

    /* creates the output array */
    out = fftw_malloc(sizeof(fftw_complex) * size);

    /* let's gooooo */
    plan_forward = fftw_plan_dft(d - 1, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    /* finds the maximum */
    best = 0.0;
    idx = 0;
    for (i = 0; i < size; ++i) {
        if (out[i][0] > best) {
            best = out[i][0];
            idx = i;
        }
    }

    /* recovers the vector and puts it in v */
    for (i = 0; i < d - 1; ++i) {
        v[d - 2 - i].value = idx % q;
        idx = (idx - v[d - 2 - i].value)/q;
    }

    /* frees the memory for fft */
    free(n);
    fftw_destroy_plan(plan_forward);
    fftw_free(in);
    fftw_free(out);
}

void bkw_free_log() {
    free(log_j);
    log_j = NULL;
}
