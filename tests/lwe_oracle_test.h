/*
 * lwe_oracle_test.h : Unit tests for the LWE oracle.
 *
 *  Created on: Apr 26, 2016
 *      Author: Aymeric Genet
 */

/* nametag to call test at execution */
#define LWE_ORACLE_TEST_NAME "lwe_oracle"

/*
 * Tests lwe_oracle.lwe_oracle_predef()
 */

char * test_lwe_oracle_predef();

/*
 * Tests lwe_oracle.lwe_oracle()
 */

char * test_lwe_oracle();

/*
 * Gathers all tests in one function.
 */

char * lwe_oracle_all_tests();
