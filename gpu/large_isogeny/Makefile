CXX ?= g++
NVCC ?= nvcc
ARCH ?= sm_90

.PHONY: test gpu clean
test: test_cpu
	./test_cpu

test_cpu: test_cpu.cpp velu.hpp
	$(CXX) -O2 -std=c++17 -Wall -Wextra -Werror -o $@ $<

gpu: test_gpu
batch.o: batch.cu velu.hpp
	$(NVCC) -O3 -std=c++17 -arch=$(ARCH) -c -o $@ $<

test_gpu: test_gpu.cu batch.o velu.hpp
	$(NVCC) -O3 -std=c++17 -arch=$(ARCH) -o $@ test_gpu.cu batch.o

gpu-test: test_gpu
	./test_gpu

clean:
	rm -f test_cpu batch.o test_gpu
