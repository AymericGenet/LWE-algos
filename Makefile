
main:
	gcc -o main $(wildcard src/*) main.c -g -Wall -ansi -pedantic -lm -lfftw3

tests:
	gcc -o test $(wildcard src/*) $(wildcard tests/*) tests.c -g -Wall -ansi -pedantic -lm -lfftw3

bench:
	gcc -o benchmark $(wildcard src/*) $(wildcard bench/*) bench.c -g -Wall -ansi -pedantic -lm -lfftw3

scratchpad:
	gcc -o output $(wildcard src/*) scratchpad.c -g -Wall -ansi -pedantic -lm -lfftw3

.PHONY: tests scratchpad main bench
