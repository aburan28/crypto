# gpu/fes -- build + test the Gray-code FES search kernels.
#
#   make test            host verification (no GPU): the same fes.cuh the
#                        CUDA/Metal kernels use, checked against brute force
#   make test N=20 M=22 SEED=3   a different generated system
#   make cuda            compile+link the CUDA launcher (needs nvcc)
#   make nvcc-check      compile-only for several SM architectures
#   make metal           compile the Metal shader (needs macOS + xcrun)
#   make clean
#
# `make test` on a GPU-less machine proves the algorithm: the field-free
# quadratic evaluation, the derivative-maintained Gray walk, and the sub-cube
# sharding that maps one thread to one sub-cube.  Only the launch is
# device-specific.

PY        ?= python3
CXX       ?= g++
NVCC      ?= nvcc
CXXFLAGS  ?= -O2 -std=c++17 -Wall -Wextra
ARCH      ?= sm_90
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo

N    ?= 18
M    ?= 20
SEED ?= 1

GEN  := system.h
HDRS := fes.cuh

.PHONY: all test cuda nvcc-check metal clean gen
all: test

gen: $(GEN)

system.h: fref.py
	$(PY) fref.py system --n $(N) --m $(M) --seed $(SEED) > $@

test_cpu_main: test_cpu.cpp $(HDRS) $(GEN)
	$(CXX) $(CXXFLAGS) -I. -o $@ $<

test: test_cpu_main
	./test_cpu_main

cuda: fes.cu $(HDRS) $(GEN)
	$(NVCC) $(NVCCFLAGS) -o fes_cuda fes.cu

# Compile-only sweep across architectures, like gpu/ecc2k's nvcc-check.
nvcc-check: fes.cu $(HDRS) $(GEN)
	@for a in sm_70 sm_80 sm_90; do \
	  echo "=== nvcc -arch=$$a (compile only) ==="; \
	  $(NVCC) -O3 -std=c++17 -arch=$$a -Xptxas -v -c fes.cu -o /dev/null || exit 1; \
	done

# Metal shader compile-check (macOS only).
metal: fes.metal
	xcrun -sdk macosx metal -std=metal3.0 -c fes.metal -o fes.air

clean:
	rm -f test_cpu_main fes_cuda fes.air system.h
