CXX ?= g++
NVCC ?= nvcc
ARCH ?= sm_120
NIST_PTX ?= 1
MAXREG ?= 0
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unknown-pragmas
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo -Xptxas=-v
ifneq ($(MAXREG),0)
NVCCFLAGS += -maxrregcount=$(MAXREG)
endif

.PHONY: all test bench ptx clean
all: test
test: test_cpu
	./test_cpu
test_cpu: test_cpu.cpp nist_fields.cuh nist_ptx.cuh nist_solinas.cuh
	$(CXX) $(CXXFLAGS) -o $@ test_cpu.cpp
bench: bench.cu nist_fields.cuh nist_ptx.cuh nist_solinas.cuh
	$(NVCC) $(NVCCFLAGS) -DNIST_PTX=$(NIST_PTX) -o $@ bench.cu
ptx: bench.cu nist_fields.cuh nist_ptx.cuh nist_solinas.cuh
	$(NVCC) $(NVCCFLAGS) -DNIST_PTX=$(NIST_PTX) -ptx -o bench.ptx bench.cu
clean:
	rm -f test_cpu bench bench.ptx
