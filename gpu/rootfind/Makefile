# gpu/rootfind -- root-finding 2-decomposition oracle: CUDA + host emulation.
#
#   make test            CPU verification of the shared headers (no GPU, no nvcc)
#   make lib             librootfind_oracle.a, emulation backend only
#                        (what gpu/rootfind/sys/build.rs builds for cargo --features gpu-rootfind)
#   make lib-cuda        librootfind_cuda.a, emulation + CUDA backend (needs nvcc)
#   make bench           CUDA driver: selftest vs emulation + throughput (needs nvcc)
#   make bench ARCH=sm_100      datacenter Blackwell; default sm_120 (RTX PRO 6000 Blackwell)
#   make ptxas-info      per-kernel registers / spills / smem for the chosen ARCH
#   make clean
#
# Conventions follow gpu/ecc: every arithmetic function is
# `__host__ __device__`, so `make test` compiles the *same* headers the
# kernels use.  Blackwell flags: -gencode arch=compute_120,code=sm_120 for
# native SASS plus -gencode arch=compute_120,code=compute_120 for a PTX
# fallback that JIT-compiles on future architectures.

CXX      ?= c++
NVCC     ?= nvcc
AR       ?= ar
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas -pthread
ARCH     ?= sm_120
# compute_XX is derived from sm_XX
COMPUTE  := $(subst sm_,compute_,$(ARCH))
GENCODE  ?= -gencode arch=$(COMPUTE),code=$(ARCH) -gencode arch=$(COMPUTE),code=$(COMPUTE)
NVCCFLAGS ?= -O3 -std=c++17 -lineinfo $(GENCODE) -Xcompiler -pthread
# Launch-bounds knobs (see rootfind_cuda.cu header): block size and the
# minimum resident blocks per SM ptxas must honour (0 = unconstrained).
RF_BLOCK      ?= 128
RF_MIN_BLOCKS ?= 0
RF_DEFS := -DRF_BLOCK=$(RF_BLOCK) -DRF_MIN_BLOCKS=$(RF_MIN_BLOCKS)

HDRS := fp64.cuh fp256x64.cuh rootfind.cuh rootfind_host.hpp rootfind_api.h

.PHONY: all test lib lib-cuda bench ptxas-info clean

all: test

test_cpu: test_cpu.cpp rootfind_api.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -DRF_HAVE_CUDA=0 -o $@ test_cpu.cpp rootfind_api.cpp

test: test_cpu
	./test_cpu

# --- emulation-only static library (no CUDA toolchain required) -------
rootfind_api_emul.o: rootfind_api.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -DRF_HAVE_CUDA=0 -c -o $@ rootfind_api.cpp

librootfind_oracle.a: rootfind_api_emul.o
	$(AR) rcs $@ $^

lib: librootfind_oracle.a

# --- CUDA-enabled static library --------------------------------------
rootfind_api_cuda.o: rootfind_api.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -DRF_HAVE_CUDA=1 -c -o $@ rootfind_api.cpp

rootfind_cuda.o: rootfind_cuda.cu $(HDRS)
	$(NVCC) $(NVCCFLAGS) $(RF_DEFS) -DRF_HAVE_CUDA=1 -c -o $@ rootfind_cuda.cu

librootfind_cuda.a: rootfind_api_cuda.o rootfind_cuda.o
	$(AR) rcs $@ $^

lib-cuda: librootfind_cuda.a

# --- CUDA driver -------------------------------------------------------
bench: bench.cu rootfind_api_cuda.o rootfind_cuda.o $(HDRS)
	$(NVCC) $(NVCCFLAGS) $(RF_DEFS) -DRF_HAVE_CUDA=1 -o $@ bench.cu rootfind_api_cuda.o rootfind_cuda.o

ptxas-info: rootfind_cuda.cu $(HDRS)
	$(NVCC) $(NVCCFLAGS) $(RF_DEFS) -DRF_HAVE_CUDA=1 -Xptxas -v -c -o /dev/null rootfind_cuda.cu

clean:
	rm -f test_cpu bench *.o *.a
