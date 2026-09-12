# gpu/macaulay -- build + test.
#
#   make test          generated headers + CPU verification (no GPU needed)
#   make bench         CUDA benchmark binary (needs nvcc)
#   make bench ARCH=sm_100     datacenter Blackwell; sm_120 for RTX 50xx
#   make clean
#
# The CPU tests compile the *same* .cuh headers with g++, so a green
# `make test` on a GPU-less machine means the field arithmetic and the
# serial reduction are right.  The block-parallel kernel is device-only
# and is checked by `./bench selftest` against that verified serial path.

PY       ?= python3
CXX      ?= g++
NVCC     ?= nvcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter
ARCH     ?= sm_90
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo

# p and the shape default to the measured gaudry_cubic case: the largest
# size RESEARCH_RESIDUAL_WALKS.md section 11.6 runs, and the 226 x 286
# matrix its row selection produces.
P      ?= 2083
ROWS   ?= 226
COLS   ?= 286
SEED   ?= 0x9E3779B9
DEFICIT?= 4

GEN   := params.h vectors.h
HDRS  := fpmod.cuh macaulay.cuh
TESTS := test_cpu_main

.PHONY: all test bench clean gen
all: test

gen: $(GEN)

params.h: mref.py
	$(PY) mref.py params --p $(P) --rows $(ROWS) --cols $(COLS) > $@

vectors.h: mref.py params.h
	$(PY) mref.py vectors --p $(P) --rows $(ROWS) --cols $(COLS) \
	  --seed $(SEED) --deficit $(DEFICIT) > $@

test_cpu_main: test_cpu.cpp $(HDRS) $(GEN)
	$(CXX) $(CXXFLAGS) -I. -o $@ $<

test: $(TESTS)
	@rc=0; for t in $(TESTS); do echo "=== ./$$t ==="; ./$$t || rc=1; done; exit $$rc

# A second shape, to catch anything that only works at the default one.
.PHONY: test-small
test-small:
	$(MAKE) clean
	$(MAKE) test P=271 ROWS=64 COLS=80
	$(MAKE) clean
	$(MAKE) test

bench: bench.cu $(HDRS) $(GEN)
	$(NVCC) $(NVCCFLAGS) -I. -o $@ $<

clean:
	rm -f $(GEN) $(TESTS) bench
