/*
 * bkw_test.h : Unit tests for the BKW algorithm.
 *
 *  Created on: Apr 28, 2016
 *      Author: Aymeric Genet
 */

/* nametag to call test at execution */
#define BKW_TEST_NAME "bkw"


/*
 * Tests bkw.bkw_lf1()
 */

char * test_bkw_lf1();

/*
 * Tests bkw.bkw_lf2()
 */

char * test_bkw_lf2();

/*
 * Shows extracted noise distribution.
 */

char * test_bkw_distrib();

/*
 * Tests bkw.bkw_hypo_testing()
 */

char * test_bkw_hypo_testing();

/*
 * Tests bkw.bkw_fft()
 */

char * test_bkw_fft();

/*
 *  Gathers all tests in one function.
 */

char * bkw_all_tests();
