VPATH := cpp fortran

.PHONY:
all: bin/poisson-nvcpp  \
     bin/poisson-nvfort-serial \
     bin/poisson-nvfort-multicore \
     bin/poisson-icpc \
     bin/poisson-ifort-serial \
     bin/poisson-ifort-multicore

.PHONY:
clean:
	rm -rf bin build

bin/poisson-nvcpp: poisson.cpp | bin
	nvc++ -O3 -std=c++20 -mp $^ -o $@

bin/poisson-nvfort: utils.f90 poisson.f90 | bin
	nvfortran -O3 -stdpar=multicore $^ -o $@

bin/poisson-nvfort-multicore: utils.f90 poisson.f90 | bin
	nvfortran -O3 -stdpar=multicore $^ -o $@

bin/poisson-icpc: poisson.cpp | bin
	icpc -O3 -std=c++20 -qopenmp $^ -o $@

bin/poisson-ifort: utils.f90 poisson.f90 | bin
	ifort -module build -O3 $^ -o $@

bin/poisson-ifort-multicore: utils.f90 poisson.f90 | bin
	ifort -module build -O3 -parallel $^ -o $@

bin:
	mkdir -p bin build
