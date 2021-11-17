poisson: poisson.cpp
	nvc++ -O3 -std=c++20 $^ -o $@
