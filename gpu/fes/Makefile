# gpu/fes -- build + test the Gray-code FES search kernels.
#
#   make test            host verification (no GPU): the same fes.cuh the
#                        CUDA/Metal kernels use, checked against brute force
#   make test N=20 M=22 SEED=3   a different generated system
#   make host-worker     the g++ FES worker (fes_solve) the Rust binary drives
#   make cuda            compile+link the CUDA worker (needs nvcc)
#   make nvcc-check      compile-only for several SM architectures
#   make metal           compile the Metal shader (needs macOS + xcrun)
#   make metal-worker    the Metal FES worker binary (needs macOS + clang++)
#   make clean
#
# fes_solve, fes_cuda and fes_metal all speak the same file I/O contract
# (fes_io.hpp), so the Rust `icx` binary drives whichever is present through
# one code path and re-verifies every returned solution.
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
HDRS := fes.cuh fes_io.hpp

.PHONY: all test host-worker cuda nvcc-check metal metal-worker clean gen
all: test host-worker

gen: $(GEN)

system.h: fref.py
	$(PY) fref.py system --n $(N) --m $(M) --seed $(SEED) > $@

test_cpu_main: test_cpu.cpp fes.cuh $(GEN)
	$(CXX) $(CXXFLAGS) -I. -o $@ $<

test: test_cpu_main
	./test_cpu_main

# The host-emulation FES worker the Rust binary drives (and the CPU fallback).
# No GPU and no generated system.h needed — it reads the system at run time.
host-worker: fes_solve
fes_solve: fes_solve.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -I. -o $@ $<

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

# The Metal FES worker binary (macOS only): Objective-C++ host + the shader.
metal-worker: fes_metal.mm fes.metal fes_io.hpp
	xcrun -sdk macosx metal -std=metal3.0 -c fes.metal -o fes.air
	xcrun -sdk macosx metallib fes.air -o fes.metallib
	clang++ $(CXXFLAGS) -I. -ObjC++ -framework Metal -framework Foundation \
	  -o fes_metal fes_metal.mm

clean:
	rm -f test_cpu_main fes_solve fes_cuda fes_metal fes.air fes.metallib system.h
