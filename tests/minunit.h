/*
 * minunit.h : Minimal unit testing
 *
 * Source: http://www.jera.com/techinfo/jtns/jtn002.html
 */

#include <stdio.h>

extern int tests_run;

#define mu_assert(message, test) do { \
    if (!(test)) { \
        return message; \
    } \
} while (0)

#define mu_run_test(test, name) do { \
    char * message; \
    test_head(name); \
    message = test(); \
    tests_run++; \
    if (message) { \
        return message; \
    } \
    test_end; \
} while (0)

#define test_head(name) printf("Test %s:\n",name);
#define test_end printf("\t[OK]\n");
