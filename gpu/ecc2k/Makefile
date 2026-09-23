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
GEN    := $(foreach c,$(CURVES),curve_$(c).h vec_$(c).h) vec_canon.h vec_fold.h
TESTS  := test_ecc2k95 test_k23 test_k41 \
          $(foreach c,$(CURVES),test_pt_$(c)) \
          $(foreach c,$(CURVES),test_emu_$(c))

HDRS := f2m.cuh koblitz.cuh rho2k.cuh scalar_ring.hpp rho2k_host.hpp \
        pairtable.cuh

.PHONY: all test bench clean gen roundtrip nvcc-check
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

# A folded pair table as `PairSumTable::build_folded_within` stored it, and
# the row plan it was built from, so `pairtable_fold_kernel` can be checked
# against the table it has to reproduce rather than against a second
# account of what that table is.  Explicit for the same reason as above.
vec_fold.h: ../../examples/dump_fold_vectors.rs ../../src/cryptanalysis/koblitz_fast.rs \
            ../../src/cryptanalysis/koblitz_index_calculus.rs
	cd ../.. && $(CARGO) run --release --quiet --example dump_fold_vectors \
	  > gpu/ecc2k/vec_fold.h

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

# `pairtable_kernel` and `pairtable_fold_kernel` themselves, compiled as
# host C++ behind a `__CUDACC__` shim and run once per emulated thread.
# Everything above tests what the kernels are made of; this is the only
# thing that tests the kernels, and the fold against the CPU's own table.
test_emu_$(1): test_pairtable_emu.cpp $$(HDRS) curve_$(1).h vec_canon.h vec_fold.h
	$$(CXX) $$(CXXFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_$(1).h"' -o $$@ $$<
endef
$(foreach c,$(CURVES),$(eval $(call TEST_RULE,$(c))))

test: $(TESTS)
	@rc=0; for t in $(TESTS); do \
	  echo "=== ./$$t ==="; ./$$t || rc=1; done; \
	echo "=== roundtrip ==="; $(MAKE) --no-print-directory roundtrip || rc=1; \
	exit $$rc

# The fold's round trip: the table the emulated kernel built, written out
# and taken over by `PairSumTable::from_folded_parts`, then asked every
# question the CPU's own table answers.  `test_emu_*` compares stored
# state; this compares answers, and is what runs every tag through the
# recovery that reads it.
FOLD_TABLES := fold_table_k23.bin fold_table_k41.bin
fold_table_%.bin: test_emu_%
	./test_emu_$* $@ > /dev/null

roundtrip: $(FOLD_TABLES)
	cd ../.. && $(CARGO) run --release --quiet --example load_fold_table -- \
	  $(addprefix gpu/ecc2k/,$(FOLD_TABLES))

# Compile, don't run: every kernel here through nvcc.  Nothing else does --
# the host tests build `pairtable.cuh`'s kernels with g++ behind a
# `__CUDACC__` shim, and `ptx_stats2k.sh` compiles the walk with clang --
# so this is the one check that the code a device would run gets past the
# compiler a device build would use.  No GPU needed.  Warnings are errors,
# and ptxas reports what each pair-table kernel takes in registers and
# stack.  `bench2k.cu` stands in for `kernels2k.cuh`, which it includes.
NVCC_ARCHS ?= sm_90 sm_100 sm_120
NVCC_CHECK_FLAGS := -O3 -std=c++17 -Werror all-warnings
nvcc-check: $(HDRS) kernels2k.cuh bench2k.cu curve_k23.h curve_k41.h curve_ecc2k95.h
	@for a in $(NVCC_ARCHS); do \
	  for c in k23 k41; do \
	    echo "=== nvcc -arch=$$a pairtable.cuh, $$c ==="; \
	    out=$$($(NVCC) $(NVCC_CHECK_FLAGS) -arch=$$a -Xptxas -v -x cu -c -o /dev/null \
	      -DGPU_ECC2K_CURVE_HEADER="\"curve_$$c.h\"" pairtable.cuh 2>&1) \
	      || { echo "$$out"; exit 1; }; \
	    echo "$$out" | grep -E "entry function|Used"; \
	  done; \
	  echo "=== nvcc -arch=$$a bench2k.cu, ecc2k95 ==="; \
	  $(NVCC) $(NVCC_CHECK_FLAGS) -arch=$$a -c -o /dev/null \
	    -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' bench2k.cu || exit 1; \
	done

bench: bench2k.cu $(HDRS) kernels2k.cuh curve_ecc2k95.h
	$(NVCC) $(NVCCFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' -o $@ $<

clean:
	rm -f $(GEN) $(TESTS) $(FOLD_TABLES) bench
