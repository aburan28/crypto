# gpu/btcpuzzle -- Pollard kangaroo for interval ECDLP on secp256k1.
#
#   make test      CPU verification (no GPU needed)
#   make bench     CUDA solver + benchmark (needs nvcc)
#   make clean
#
# Reuses the verified secp256k1 field and point arithmetic from ../ecc.

PY       ?= python3
CXX      ?= g++
NVCC     ?= nvcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas -I../ecc
ARCH     ?= sm_90
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo -I../ecc

HDRS := kangaroo.cuh kangaroo_host.hpp ../ecc/fp256.cuh ../ecc/point.cuh

.PHONY: all test bench clean
all: test

../ecc/curve_secp256k1.h:
	$(MAKE) -C ../ecc curve_secp256k1.h

test_kangaroo: test_kangaroo.cpp $(HDRS) ../ecc/curve_secp256k1.h
	$(CXX) $(CXXFLAGS) -o $@ $<

test: test_kangaroo
	./test_kangaroo

# Slow: 120 solves per configuration.  Not part of `make test`.
tune_kangaroo: tune_kangaroo.cpp $(HDRS) ../ecc/curve_secp256k1.h
	$(CXX) $(CXXFLAGS) -o $@ $<

bench: kangaroo.cu kernels_kangaroo.cuh $(HDRS) ../ecc/curve_secp256k1.h
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -f test_kangaroo tune_kangaroo bench
