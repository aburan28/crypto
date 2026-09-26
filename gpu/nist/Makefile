CXX ?= g++
NVCC ?= nvcc
ARCH ?= sm_120
NIST_PTX ?= 1
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo -Xptxas=-v

.PHONY: all test bench ptx clean
all: test
test: test_cpu
	./test_cpu
test_cpu: test_cpu.cpp nist_fields.cuh
	$(CXX) $(CXXFLAGS) -o $@ test_cpu.cpp
bench: bench.cu nist_fields.cuh nist_ptx.cuh
	$(NVCC) $(NVCCFLAGS) -DNIST_PTX=$(NIST_PTX) -o $@ bench.cu
ptx: bench.cu nist_fields.cuh nist_ptx.cuh
	$(NVCC) -O3 -std=c++17 -arch=$(ARCH) -DNIST_PTX=$(NIST_PTX) -ptx -o bench.ptx bench.cu
clean:
	rm -f test_cpu bench bench.ptx
