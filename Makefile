
main:
	gcc -o main $(wildcard src/*) main.c -g -Wall -ansi -pedantic -lm

tests:
	gcc -o test $(wildcard src/*) $(wildcard tests/*) tests.c -g -Wall -ansi -pedantic -lm

scratchpad:
	gcc -o output $(wildcard src/*) scratchpad.c -g -Wall -ansi -pedantic -lm

.PHONY: tests scratchpad main
