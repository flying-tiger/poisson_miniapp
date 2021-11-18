all: poisson-nvcpp poisson-icpc

poisson-nvcpp: poisson.cpp
	nvc++ -O3 -std=c++20 $^ -o $@

poisson-icpc: poisson.cpp
	icpc -O3 -std=c++20 $^ -o $@
