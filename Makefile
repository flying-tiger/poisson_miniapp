.PHONY:
all: poisson-nvcpp poisson-nvfort poisson-icpc poisson-ifort

.PHONY:
clean:
	rm -f *.o *.mod poisson-*

poisson-nvcpp: poisson.cpp
	nvc++ -O3 -std=c++20 $^ -o $@

poisson-nvfort: utils.f90 poisson.f90
	nvfortran -O3 $^ -o $@

poisson-icpc: poisson.cpp
	icpc -O3 -std=c++20 $^ -o $@

poisson-ifort: utils.f90 poisson.f90
	ifort -O3 $^ -o $@
