Introduction
============

The Learning With Errors (LWE) problem is a fundamental problem cryptography as it is used in the construction of many cryptographic primitives (e.g. in fully homomorphic encryption). In the recent years there have been several papers that study the complexity of solving LWE.

Albrecht et al. adapted the BKW algorithm to solve LWE and introduced the lazy modulus switching. Later, Duc et al. improved these results by using multidimensional Fourier Transform.

The goal of this project is to implement the existing algorithms and see the difference between the theoretical results and the ones obtained in practice.

Configuration
=============

The project works with the following versions of the tools :

* gcc (GCC) 4.8.1
* GNU Make 3.81
* FFTW 3.3.4

How to launch
=============

To build the main program, type

    $> make main

Then launch it with the following parameters :

    -n size       size for the secret s (default: 6)
    -q prime      modulus for ring Z_q (default: 37)
    -a depth      oracle depth (default: 5)
    -m samples    gives an amount of samples (default: 15)

    -i index      launches with precomputed parameters chosen by index (default: 0)

    -r range      launches over specified range (default: 1000)

    -x noise      chooses noise distribution (default: 0)
                    0=Rounded Gaussian noise
                    1=Discrete Gaussian noise
                    2=Discrete Uniform noise
    -l algo       chooses sample reduction algorithm (default: 0)
                    0=LF1
                    1=LF2
    -t algo       chooses solving algorithm (default: 0)
                    0=Log-likelihood
                    1=Fast Fourier tranfsorm

    -v            verboses wrong guesses along with extracted noise distribution

    -h            prints this



References
==========

In this project, we make use of the following source code :

* JTN002 - MinUnit -- a minimal unit testing framework for C : http://www.jera.com/techinfo/jtns/jtn002.html
* benchmarking - How can I benchmark C code easily? - Stack Overflow : http://stackoverflow.com/questions/2349776/how-can-i-benchmark-c-code-easily
* Standalone erf function : http://www.johndcook.com/blog/cpp_erf/
* bkw-lwe : https://bitbucket.org/malb/bkw-lwe/src

We implemented concepts, ideas, and algorithms that come from the following papers :

* Albrecht, M. et al. 2015. *On the complexity of the BKW algorithm*
* Duc, A. et al.      2015. *Better algorithms for LWE and LWR*
* Choi, G. et al.     2016. *Study LWE algorithms*



Future work
===========

 * Lazy Modulus switching



Contact
=======

Author : Aymeric Genet <aymeric {dot} genet {at} epfl {dot} ch>
Supervisor: Sonia Mihaela Bogos <soniamihaela {dot} bogos {at} epfl {dot} ch>
Professor:  Serge Vaudenay <serge {dot} vaudenay {at} epfl {dot} ch>



Thanks for reading !
