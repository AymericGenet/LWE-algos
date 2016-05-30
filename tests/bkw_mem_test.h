/*
 * bkw_mem_test.h
 *
 *  Created on: May 29, 2016
 *      Author: Aymeric Genet
 */

/* nametag to call test at execution */
#define BKW_MEM_TEST_NAME "bkw_mem"


char * test_bkw_mem_lf1();

char * test_bkw_mem_lf2();

char * test_bkw_mem_hypo_testing();

char * test_bkw_mem_fft();

/*
 *  Gathers all tests in one function.
 */

char * bkw_mem_all_tests();
