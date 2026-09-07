# gpu/ecc2k -- build + test.
#
#   make test          generated headers + CPU verification (no GPU needed)
#   make bench         CUDA benchmark binary (needs nvcc)
#   make bench ARCH=sm_100
#   make clean

PY      ?= python3
CXX     ?= g++
NVCC    ?= nvcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas
ARCH    ?= sm_90
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo

CURVES := ecc2k95 k23 k41
GEN    := $(foreach c,$(CURVES),curve_$(c).h vec_$(c).h)
TESTS  := test_ecc2k95 test_k23 test_k41

HDRS := f2m.cuh koblitz.cuh rho2k.cuh scalar_ring.hpp rho2k_host.hpp

.PHONY: all test bench clean gen
all: test

gen: $(GEN)

curve_%.h: ecref2k.py
	$(PY) ecref2k.py params --curve $* > $@
vec_%.h: ecref2k.py
	$(PY) ecref2k.py vectors --curve $* > $@

define TEST_RULE
test_$(1): test_cpu2k.cpp $$(HDRS) curve_$(1).h vec_$(1).h
	$$(CXX) $$(CXXFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_$(1).h"' \
	  -DGPU_ECC2K_VECTORS_HEADER='"vec_$(1).h"' -o $$@ $$<
endef
$(foreach c,$(CURVES),$(eval $(call TEST_RULE,$(c))))

test: $(TESTS)
	@rc=0; for t in $(TESTS); do \
	  echo "=== ./$$t ==="; ./$$t || rc=1; done; exit $$rc

bench: bench2k.cu $(HDRS) kernels2k.cuh curve_ecc2k95.h
	$(NVCC) $(NVCCFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' -o $@ $<

clean:
	rm -f $(GEN) $(TESTS) bench
