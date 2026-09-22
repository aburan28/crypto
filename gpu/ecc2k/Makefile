# gpu/ecc2k -- build + test.
#
#   make test          generated headers + CPU verification (no GPU needed)
#   make bench         CUDA benchmark binary (needs nvcc)
#   make bench ARCH=sm_100
#   make clean

PY      ?= python3
CARGO   ?= cargo
CXX     ?= g++
NVCC    ?= nvcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas
ARCH    ?= sm_90
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -lineinfo

CURVES := ecc2k95 k23 k41
GEN    := $(foreach c,$(CURVES),curve_$(c).h vec_$(c).h) vec_canon.h
TESTS  := test_ecc2k95 test_k23 test_k41 \
          $(foreach c,$(CURVES),test_pt_$(c)) \
          $(foreach c,$(CURVES),test_emu_$(c))

HDRS := f2m.cuh koblitz.cuh rho2k.cuh scalar_ring.hpp rho2k_host.hpp \
        pairtable.cuh

.PHONY: all test bench clean gen
all: test

gen: $(GEN)

curve_%.h: ecref2k.py
	$(PY) ecref2k.py params --curve $* > $@
vec_%.h: ecref2k.py
	$(PY) ecref2k.py vectors --curve $* > $@

# The Frobenius canon's basis change, and the keys `FrobeniusCanon` makes
# from it, so `pt_canon` can be checked against the CPU it has to agree
# with.  Generated from the Rust that owns the basis rather than built
# here: the normal element comes from a randomised search, so a second
# search would find a different basis and name the same orbits
# incompatibly.  An explicit rule, because the `vec_%.h` pattern above
# would otherwise ask `ecref2k.py` for a curve called "canon".
vec_canon.h: ../../examples/dump_canon_vectors.rs ../../src/cryptanalysis/koblitz_fast.rs
	cd ../.. && $(CARGO) run --release --quiet --example dump_canon_vectors \
	  > gpu/ecc2k/vec_canon.h

define TEST_RULE
# F2M_COUNT_OPS turns on the field-operation counters, so the per-step cost
# table in README.md is counted by the harness rather than derived by hand.
test_$(1): test_cpu2k.cpp $$(HDRS) curve_$(1).h vec_$(1).h
	$$(CXX) $$(CXXFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_$(1).h"' \
	  -DGPU_ECC2K_VECTORS_HEADER='"vec_$(1).h"' -DF2M_COUNT_OPS=1 -o $$@ $$<

# The |F|^2 pair table: batched row sums against unbatched addition.
test_pt_$(1): test_pairtable.cpp $$(HDRS) curve_$(1).h vec_$(1).h vec_canon.h
	$$(CXX) $$(CXXFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_$(1).h"' \
	  -DGPU_ECC2K_VECTORS_HEADER='"vec_$(1).h"' -o $$@ $$<

# `pairtable_kernel` itself, compiled as host C++ behind a `__CUDACC__`
# shim and run once per emulated thread.  Everything above tests what the
# kernel is made of; this is the only thing that tests the kernel.
test_emu_$(1): test_pairtable_emu.cpp $$(HDRS) curve_$(1).h vec_canon.h
	$$(CXX) $$(CXXFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_$(1).h"' -o $$@ $$<
endef
$(foreach c,$(CURVES),$(eval $(call TEST_RULE,$(c))))

test: $(TESTS)
	@rc=0; for t in $(TESTS); do \
	  echo "=== ./$$t ==="; ./$$t || rc=1; done; exit $$rc

bench: bench2k.cu $(HDRS) kernels2k.cuh curve_ecc2k95.h
	$(NVCC) $(NVCCFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' -o $@ $<

clean:
	rm -f $(GEN) $(TESTS) bench
