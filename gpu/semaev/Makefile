# gpu/semaev -- build + test.
#
#   make test           generated headers + CPU verification (no GPU)
#   make test L=4       a different subspace dimension
#   make ladder         l = 3,4,5 in sequence
#   make bench2         CUDA binary (needs nvcc)
#   make clean
#
# The CPU tests compile the *same* .cuh headers with g++, so a green
# `make test` on a GPU-less machine means the field arithmetic, the
# subspace polynomial, the quartic and the sweep are right.  Only the
# launch structure is device-specific, and `./bench2 selftest` covers it.

PY       ?= python3
CXX      ?= g++
NVCC     ?= nvcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter
ARCH     ?= sm_90
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo

# n = 3l keeps the decomposition probability near 1/3! as the factor
# base grows, which is the ladder RESEARCH_SEMAEV_DECOMPOSITION.md runs.
L   ?= 5
YES ?= 4
NO  ?= 4

GEN   := params.h vectors.h
HDRS  := gf2n.cuh decomp.cuh
TESTS := test_cpu_main

.PHONY: all test ladder bench2 clean gen
all: test

gen: $(GEN)

params.h: sref.py
	$(PY) sref.py params --l $(L) > $@

vectors.h: sref.py params.h
	$(PY) sref.py vectors --l $(L) --yes $(YES) --no $(NO) > $@

test_cpu_main: test_cpu.cpp $(HDRS) $(GEN)
	$(CXX) $(CXXFLAGS) -I. -o $@ $<

test: $(TESTS)
	@rc=0; for t in $(TESTS); do echo "=== ./$$t ==="; ./$$t || rc=1; done; exit $$rc

# Several dimensions, because the subspace polynomial's degree and the
# reduction's shape both depend on l and n.
ladder:
	@for l in 3 4 5; do \
	  $(MAKE) --no-print-directory clean >/dev/null; \
	  $(MAKE) --no-print-directory test L=$$l || exit 1; \
	done
	@$(MAKE) --no-print-directory clean >/dev/null

bench2: bench2.cu $(HDRS) $(GEN)
	$(NVCC) $(NVCCFLAGS) -I. -o $@ $<

clean:
	rm -f $(GEN) $(TESTS) bench2
