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
#include <float.h>


static double * log_j = NULL;
unsigned long lwe_oracle_calls = 0;

void bkw_create(bkw_t * bkw, lwe_t lwe, int a, int b, int d, long m) {
    size_t i;

    bkw->lwe = lwe;
    bkw->a = a;
    bkw->b = b;
    bkw->d = d;
    bkw->m = m;

    bkw->tab.first = malloc(sizeof(node_t));
    bkw->tab.states = malloc(a * sizeof(int));
    bkw->tab.sample = malloc(a * sizeof(vec_t));

    for (i = 0; i < a; ++i) {
        bkw->tab.states[i] = 0;
        bkw->tab.sample[i] = NULL;
    }

    bkw_create_node(bkw->tab.first, a, lwe.q, bkw->b);
}

void bkw_free(bkw_t * bkw) {
    size_t i;

    for (i = 0; i < bkw->a; ++i) {
        free(bkw->tab.sample[i]);
    }

    free(bkw->tab.sample);
    free(bkw->tab.states);
    bkw_free_node(bkw->tab.first, bkw->a, bkw->lwe.q, bkw->b);
}

void bkw_create_node(node_t * node, int a, long q, int b) {
    size_t i, j;
    long depth;

    depth = pow(q, b);

    node->next = NULL;
    node->T = malloc(a * sizeof(vec_t *));
    for (i = 0; i < a; ++i) {
        node->T[i] = malloc(depth * sizeof(vec_t));
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

int bkw_lf1(vec_t res, bkw_t bkw, int l, math_t ** aux) {
    int n, b, d;
    long q;
    size_t i, idx;
    int start, end;
    vec_t ** T;

    if (l == 0) {
        lwe_oracle_calls++;
        return lwe_oracle(res, bkw.lwe); /* TODO : time when verbose */
    }

    /* for visibility concerns */
    n = bkw.lwe.n;
    q = bkw.lwe.q;
    b = bkw.b;
    d = bkw.d;
    T = bkw.tab.first->T;

    start = (l - 1) * b;
    end   = (l == bkw.a) ? n - d : l * b;

    while (start >= n) {
        start -= b;
    }
    if (end > n - d) {
        end = n - d;
    }
    while (1) {
        /* samples from previous layer */
        if (!bkw_lf1(aux[0], bkw, l - 1, aux + 1)) {
            return 0; /* false */
        }

        /* if first elements already 0, returns it */
        if (zero(aux[0], start, end)) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = aux[0][i];
            }
            return 1; /* true */
        }

        /* checks if collision */
        idx = index(aux[0], q, start, end);
        if (T[l-1][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = (aux[0][i] + q - T[l-1][idx][i]) % q;
            }
            return 1; /* true */
        }

        /* negates it */
        for (i = 0; i < n + 1; ++i) {
            aux[0][i] = (q - aux[0][i]) % q;
        }

        /* checks if collision */
        idx = index(aux[0], q, start, end);
        if (T[l-1][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = (aux[0][i] + q - T[l-1][idx][i]) % q;
            }
            return 1; /* true */
        }

        /* otherwise, stores it, and repeat */
        T[l-1][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l-1][idx][i] = aux[0][i];
        }
    }
    return 0; /* should not happen */
}

int bkw_lf2(vec_t res, bkw_t bkw, int l, math_t ** aux) {
    size_t i, idx, n_idx;
    int layer, n, b; /*, d; */
    long q;
    table_t * tab;
    node_t * current;
    node_t * next;
    vec_t ** T;
    vec_t sample;

    /* if call to lwe_oracle */
    if (l == 0) {
        return lwe_oracle(res, bkw.lwe);
    }

    /* for visibility concerns */
    tab = &(bkw.tab);
    n = bkw.lwe.n;
    q = bkw.lwe.q;
    b = bkw.b;
    /* d = bkw.d; *//* TODO d reduction */

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
                res[i] = (sample[i] + q - T[l][idx][i]) % q;
            }
            return 1; /* true */
        }

        /* negates it */
        for (i = 0; i < n + 1; ++i) {
            aux[0][i] = (q - sample[i]) % q;
        }

        /* checks if collision */
        n_idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][n_idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = (aux[0][i] + q - T[l][n_idx][i]) % q;
            }
            return 1; /* true */
        }

        /* otherwise, stores it */
        T[l][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l][idx][i] = sample[i];
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
            T[l][idx][i] = sample[i];
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
        if (!bkw_lf2(aux[0], bkw, l - 1, aux + 1)) {
            return 0; /* false */
        }

        /* if first elements already 0, returns it */
        if (zero(aux[0], (l - 1) * b, l * b)) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = aux[0][i];
                sample[i] = aux[0][i];
            }
            tab->states[l] = -1;
            return 1; /* true */
        }

        /* checks if collision */
        idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = (aux[0][i] + q - T[l][idx][i]) % q;
                sample[i] = aux[0][i];
            }
            return 1; /* true */
        }

        /* negates it */
        for (i = 0; i < n + 1; ++i) {
            aux[0][i] = (q - aux[0][i]) % q;
        }

        /* checks if collision */
        n_idx = index(aux[0], q, (l - 1) * b, l * b);
        if (T[l][n_idx] != NULL) {
            for (i = 0; i < n + 1; ++i) {
                res[i] = (aux[0][i] + q - T[l][n_idx][i]) % q;
                sample[i] = aux[0][i];
            }
            return 1; /* true */
        }

        /* otherwise, stores it, and repeat */
        T[l][idx] = malloc((n + 1) * sizeof(math_t));
        for (i = 0; i < n + 1; ++i) {
            T[l][idx][i] = (q - aux[0][i]) % q;
        }
    }
    return 0;
}

void bkw_hypo_testing(vec_t v, vec_t * F, bkw_t bkw, math_t ** aux) {
    int i, d, m;
    long q;
    size_t j;
    double p, best, current;

    /* for visibility concerns */
    q = bkw.lwe.q;
    d = bkw.d + 1;
    m = bkw.m;

    /* precomputes each value of log(j) */
    if (log_j == NULL) {
        log_j = malloc(q * sizeof(double));
        for (i = 0; i < q; ++i) {
            p = distributions[bkw.lwe.distrib](i > q/2 ? i - q : i, bkw.lwe.sig, q);
            log_j[i] = log(p) - log((pow(q, d-2) - p)/(double)(pow(q, d-1) - 1.0));
            log_j[i] = log_j[i] / log(2.0);
        }
    }

    /* aux[0][i] = first coordinate of a_i */
    for (i = 0; i < m; ++i) {
        aux[0][i] = F[i][0];
    }

    best = -DBL_MAX;
    /* foreach v (aux[1]) in Z_q^d, starting at v = (1 0 0 0...) */
    for (aux[1][0] = 1; aux[1][d - 2] != q; ) {
        /* foreach sample, update S[v] with a_i - c_i */
        current = 0.0;
        for (i = 0; i < m; ++i) {
            /* S[v]_i = W(aux[0][i] - c_i) */
            current += log_j[(F[i][d - 1] - aux[0][i] + q) % q];
            /* aux[0][i] += first coordinate of a_i */
            aux[0][i] = (aux[0][i] + F[i][0]) % q;
        }
        if (current/m > best) {
            for (i = 0; i < d - 1; ++i) {
                v[i] = aux[1][i];
            }
            best = current/m;
        }

        /* increases v by (1 0 0 0...) */
        aux[1][0] += 1;

        /* when v becomes (q 0 0 0...), increases the next coordinates,
        repeats the process if next coordinates also hit q (max). */
        for (i = 0; aux[1][i] >= q && i < d - 2 ; ++i) {
            aux[1][i] = 0;
            /* if there is one coordinate after this one */
            if (i + 1 < d - 1) {
                aux[1][i + 1] += 1;
                for (j = 0; j < m; ++j) {
                    aux[0][j] = aux[0][j] + F[j][i + 1];
                }
            } else {
                aux[1][i] = q; /* to get outside of the loop */
            }
        }
    }
}

void bkw_fft(vec_t v, vec_t * F, bkw_t bkw) {
    size_t i, j, idx, size;
    long q;
    double best;
    int * n, d, m;
    fftw_complex * in;
    fftw_complex * out;
    fftw_plan plan_forward;

    /* for visibility concerns */
    q = bkw.lwe.q;
    d = bkw.d + 1;
    m = bkw.m;

    /* creates the input array */
    n = malloc((d - 1) * sizeof(int));
    size = 1;

    for (i = 0; i < d - 1; ++i) {
        size *= q;
        n[i] = q;
    }
    in = fftw_malloc(sizeof(fftw_complex) * size);

    for (i = 0; i < size; ++i) {
        in[i][0] = 0.0;
        in[i][1] = 0.0;
    }

    /* puts value to be transformed */
    for (i = 0; i < m; ++i) {
        idx = 0;
        for (j = 0; j < d - 1; ++j) {
            idx = F[i][j] + q*idx;
        }
        in[idx][0] += cos(2 * PI_VAL * F[i][d - 1] /((double) q));
        in[idx][1] += sin(2 * PI_VAL * F[i][d - 1] /((double) q));
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
        v[d - 2 - i] = idx % q;
        idx = (idx - v[d - 2 - i])/q;
    }

    /* frees the memory for fft */
    /* TODO pre-allocate... maybe? */
    free(n);
    fftw_destroy_plan(plan_forward);
    fftw_free(in);
    fftw_free(out);
}

void bkw_free_log() {
    free(log_j);
    log_j = NULL;
}
