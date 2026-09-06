# gpu/ecc -- build + test.
#
#   make test          generated headers + CPU verification (no GPU needed)
#   make bench         CUDA benchmark binary (needs nvcc)
#   make bench ARCH=sm_100   target Blackwell (default: sm_80;sm_90;sm_100)
#   make clean
#
# The CPU tests compile the *same* .cuh headers with g++, so a green
# `make test` on a GPU-less machine means the field/point/walk logic is
# correct; only the launch configuration and memory layout remain
# GPU-specific.

PY      ?= python3
CXX     ?= g++
NVCC    ?= nvcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas
# Blackwell: sm_100 (B100/B200, GB200), sm_120 (RTX 50xx / RTX PRO Blackwell).
# sm_90 = Hopper, sm_80 = Ampere.
ARCH    ?= sm_90
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo

GEN := curve_secp256k1.h curve_toy40.h \
       vec_secp256k1_mont.h vec_secp256k1_fast.h vec_toy40_mont.h

TESTS := test_secp_fast test_secp_mont test_toy_mont

.PHONY: all test bench clean gen
all: test

gen: $(GEN)

curve_%.h: ecref.py
	$(PY) ecref.py params --curve $* > $@

vec_secp256k1_mont.h: ecref.py
	$(PY) ecref.py vectors --curve secp256k1 --mode mont > $@
vec_secp256k1_fast.h: ecref.py
	$(PY) ecref.py vectors --curve secp256k1 --mode fast > $@
vec_toy40_mont.h: ecref.py
	$(PY) ecref.py vectors --curve toy40 --mode mont > $@

HDRS := fp256.cuh point.cuh rho.cuh rho_host.hpp

# secp256k1 with the special (fast) reduction -- the production path
test_secp_fast: test_cpu.cpp $(HDRS) curve_secp256k1.h vec_secp256k1_fast.h
	$(CXX) $(CXXFLAGS) -DGPU_ECC_CURVE_HEADER='"curve_secp256k1.h"' \
	  -DGPU_ECC_VECTORS_HEADER='"vec_secp256k1_fast.h"' -DFP_FAST=1 -o $@ $<

# secp256k1 forced through Montgomery -- proves the generic path on the
# same curve, and that fast/Montgomery agree on canonical values
test_secp_mont: test_cpu.cpp $(HDRS) curve_secp256k1.h vec_secp256k1_mont.h
	$(CXX) $(CXXFLAGS) -DGPU_ECC_CURVE_HEADER='"curve_secp256k1.h"' \
	  -DGPU_ECC_VECTORS_HEADER='"vec_secp256k1_mont.h"' -DFP_FAST=0 -o $@ $<

# 40-bit toy curve with a != 0: exercises the generic doubling formula and
# is small enough that the end-to-end rho solve actually finds the DLP
test_toy_mont: test_cpu.cpp $(HDRS) curve_toy40.h vec_toy40_mont.h
	$(CXX) $(CXXFLAGS) -DGPU_ECC_CURVE_HEADER='"curve_toy40.h"' \
	  -DGPU_ECC_VECTORS_HEADER='"vec_toy40_mont.h"' -DFP_FAST=0 -o $@ $<

test: $(TESTS)
	@rc=0; for t in $(TESTS); do \
	  echo "=== ./$$t ==="; ./$$t || rc=1; done; exit $$rc

bench: bench.cu $(HDRS) kernels.cuh curve_secp256k1.h
	$(NVCC) $(NVCCFLAGS) -DGPU_ECC_CURVE_HEADER='"curve_secp256k1.h"' -o $@ $<

clean:
	rm -f $(GEN) $(TESTS) bench
