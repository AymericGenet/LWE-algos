/*
 * math.c
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include "math.h"
#include "misc.h"

#include <math.h>
#include <stdlib.h>

double (*distributions[DISTRIB_NUMBER])(long, double, int, long) = {
    discrete_gaussian_pdf, rounded_gaussian_pdf, uniform_pdf };


long random_sample(distribution_t distrib, double sigma, long q) {
    double M, u;
    long sample, x;
    double (*pdf)(long x, double sigma, int n, long q);

    pdf = distributions[distrib];
    M = pdf(0, sigma, 1, q); /* we suppose this is the max value */
    do {
        read_drandom(&u);
        x = u * q;
        read_drandom(&u);
        u = u * M;
        sample = x > q/2 ? x - q : x;
    } while(u >= pdf(sample, sigma, 1, q));

    return x;
}

double discrete_gaussian_pdf(long x, double sigma, int n, long q) {
    double result, sum;
    long y;

    sigma *= sqrt(pow(2, n-2));
    result = exp(-(x*x)/(2*sigma*sigma));
    sum = 0;
    for (y = -q/2; y <= q/2; ++y) {
        sum += exp(-(y*y)/(2*sigma*sigma));
    }

    return result/sum;
}

double rounded_gaussian_cdf(double x, double sigma) {
    return 0.5 * (1 + custom_erf(x/sqrt(2*sigma*sigma)));
}

double rounded_gaussian_pdf(long x, double sigma, int n, long q) {
    double result;
    int i;

    if (n > 2) {
        sigma *= sqrt(pow(2, n-2));
    }
    result = 0.0;
    for (i = -GAUSSIAN_BOUND; i <= GAUSSIAN_BOUND; ++i) {
        result += rounded_gaussian_cdf(q*i + x + 0.5, sigma);
        result -= rounded_gaussian_cdf(q*i + x - 0.5, sigma);
    }

    return result;
}

double uniform_pdf(long x, double sigma, int n, long q) {
    double pdf = 0.0, p1, p2;
    int beta;
    size_t i;

    if (n == 0) {
        beta = (sqrt(12*sigma*sigma + 1) - 1)/2;

        return (x >= -beta && x <= beta) ? 1.0/(2*beta + 1) : 0.0;
    }

    for (i = 0; i < q; i++) {
        p1 = uniform_pdf(x + i > q/2 ? q - (x + i) : x + i, sigma, q, n - 1);
        p2 = uniform_pdf(i > q/2 ? q - i : i, sigma, q, n - 1);
        pdf += p1 * p2;
    }

    return pdf;
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
