CXX ?= c++
NVCC ?= nvcc
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pedantic

.PHONY: all host-worker cuda-worker clean

all: host-worker

host-worker: spmv_solve

cuda-worker: spmv_cuda

spmv_solve: spmv_solve.cpp spmv_io.hpp
	$(CXX) $(CXXFLAGS) -o $@ spmv_solve.cpp

spmv_cuda: spmv.cu spmv_io.hpp
	$(NVCC) -O3 -std=c++17 -o $@ spmv.cu

clean:
	rm -f spmv_solve spmv_cuda
