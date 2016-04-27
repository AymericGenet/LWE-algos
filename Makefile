
tests:
	gcc -o tests $(wildcard src/*) $(wildcard tests/*) tests.c -g -Wall -ansi -pedantic

scratchpad:
	gcc -o output $(wildcard maths/*) scratchpad.c -g -Wall -ansi -pedantic

.PHONY: tests scratchpad
