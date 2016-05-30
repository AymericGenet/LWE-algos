/*
 * tests.c : Fancy minimalist framework for testing structures.
 *
 * Source: http://www.jera.com/techinfo/jtns/jtn002.html
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tests/minunit.h"
#include "tests/lwe_test.h"
#include "tests/bkw_test.h"
#include "tests/bkw_mem_test.h"
#include "src/math.h"

#define SUITES_NUMBER 3
#define MAX_NAMETAG_LEN 128

/* array that regroups test functions to execute */
char * (*tests[SUITES_NUMBER])() = { lwe_all_tests, bkw_all_tests,
    bkw_mem_all_tests };

/* tags for specific execution */
char * tests_nametag[SUITES_NUMBER] = { LWE_ORACLE_TEST_NAME, BKW_TEST_NAME,
    BKW_MEM_TEST_NAME };

int tests_run = 0;

/*
 *  Runs all_tests() functions.
 */

int main(int argc, char *argv[]) {
    size_t i, j;
    int * indices;
    int tests_number;
    char * result;

    init_random();

    tests_number = SUITES_NUMBER;

    if (argc >= 2) {
        /* parses the test names to execute, and preserves the calling order */
        tests_number = argc - 1;
        indices = malloc(tests_number * sizeof(int));

        for (i = 0; i < tests_number; ++i) {
            indices[i] = -1;    /* not found */

            for (j = 0; indices[i] == -1 && j < SUITES_NUMBER; ++j) {
                if (strncmp(argv[i + 1], tests_nametag[j], MAX_NAMETAG_LEN) == 0) {
                    indices[i] = j;
                }
            } /* for all tests_nametag[j] */

            if (indices[i] == -1) {
                fprintf(stderr, "[ERROR] Incorrect test nametag: %s ", argv[i + 1]);
                free(indices);
                return 2;
            } /* if index still not found */

        } /* for all argv[i] */
    } else {
        /* puts all tests to be run */
        indices = malloc(tests_number * sizeof(int));

        for (i = 0; i < tests_number; ++i) {
            indices[i] = i;
        }
    }

    /* runs the tests, halts if one failure occurs */
    result = NULL;
    for (i = 0; result == NULL && i < tests_number; ++i) {
        result = (*tests[indices[i]])();
    }

    /* prints failure if previous test failed */
    if (result != NULL) {
        printf("%s\n", result);
    } else {
        printf("ALL TESTS PASSED\n");
    }
    printf("Tests run: %i\n", tests_run);

    free(indices);

    close_random();

    return result != NULL;
}
