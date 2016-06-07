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
#include <fcntl.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>

static int devRandom;

int init_random() {
    devRandom = open("/dev/urandom", O_RDONLY);
    return devRandom;
}

int read_random(math_t * dest) {
    int success;
    success = read(devRandom, dest, sizeof *dest);
    return success;
}

int read_drandom(double * dest) {
    unsigned long random;
    int success;

    success = read(devRandom, &random, sizeof(unsigned long));
    if (success < 0) {
        return success;
    }

    *dest = random/((double) ULONG_MAX); /* FIXME ? */

    return success;
}

int close_random() {
    return close(devRandom);
}

long discrete_gaussian(double sigma, long q) {
    double M, u;
    math_t x;
    long sample;

    M = discrete_gaussian_pdf(0, sigma, q);
    do {
        read_random(&x);
        read_drandom(&u);
        u = u * M;
        sample = q/2 - ((signed long) (x % q));
    } while(u >= discrete_gaussian_pdf(sample, sigma, q));
    return sample;
}

double discrete_gaussian_pdf(long x, double sigma, long q) {
    double result, sum;
    long y;

    result = exp(-(x*x)/(2*sigma*sigma));
    sum = 0;
    for (y = -q/2; y <= q/2; ++y) {
        sum += exp(-(y*y)/(2*sigma*sigma));
    }

    return result/sum;
}

long rounded_gaussian(double sigma, long q) {
    double u1, u2, z0;
    long result;

    read_drandom(&u1);
    read_drandom(&u2);
    z0 = sigma * sqrt(-2 * log(u1)) * cos(2 * PI_VAL * u2);
    result = ((int) (floor(z0 + 0.5))) % q;

    return (result < 0) ? result + q : result;
}

double rounded_gaussian_cdf(double x, double sigma) {
    return 0.5 * (1 + custom_erf(x/sqrt(2*sigma*sigma)));
}

double rounded_gaussian_pdf(long x, double sigma, long q, long k) {
    double result;
    int i;

    result = 0.0;
    for (i = -k; i <= k; ++i) {
        result += rounded_gaussian_cdf(q*i + x + 0.5, sigma);
        result -= rounded_gaussian_cdf(q*i + x - 0.5, sigma);
    }

    return result;
}

long uniform(double sigma, long q) {
    int beta;
    math_t tmp;

    beta = (sqrt(12*sigma*sigma + 1) - 1)/2;

    read_random(&tmp);

    return (((signed long) tmp) % beta) % q;
}

double uniform_pdf(long x, double sigma, long q) {
    int beta;

    beta = (sqrt(12*sigma*sigma + 1) - 1)/2;

    return (x > -beta && x < beta) ? 1/beta : 0;
}

size_t index(vec_t elem, long q, int a, int b) {
    size_t i, idx;
    unsigned long basis;

    idx = 0;
    basis = 1;
    for (i = a; i < b; ++i) {
        idx += basis*elem[i];
        basis *= q;
    }

    return idx;
}

void unindex(vec_t res, size_t idx, long q, int a, int b) {
    int i, index;
    long basis, count;

    index = idx;
    basis = pow(q, b - a - 1);
    for (i = b - a - 1; i >= 0; --i) {
        count = 0;
        while (index - count*basis > 0) {
            count++;
        }
        if (index - count*basis != 0) {
            count--;
        }
        res[i] = count;
        index -= basis*count;
        basis /= q;
    }
}

int zero(vec_t elem, int a, int b) {
    size_t i;
    for (i = a; i < b; ++i) {
        if (elem[i] != 0) {
            return 0; /* false */
        }
    }
    return 1; /* true */
}

int equals(vec_t u, vec_t v, int a, int b) {
    size_t i;
    for (i = a; i < b; ++i) {
        if (u[i] != v[i]) {
            return 0; /* false */
        }
    }
    return 1; /* true */
}

double custom_erf(double x) {
    double a1, a2, a3, a4, a5, p, t, y;
    int sign;

    a1 =  0.254829592;
    a2 = -0.284496736;
    a3 =  1.421413741;
    a4 = -1.453152027;
    a5 =  1.061405429;
    p  =  0.3275911;

    sign = 1;
    if (x < 0) {
        sign = -1;
    }
    x = x >= 0.0 ? x : -x;

    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * exp(-x*x);

    return sign * y;
}
