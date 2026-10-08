# gpu/ecc2k -- build + test.
#
#   make test          generated headers + CPU verification (no GPU needed)
#   make nvcc-check    every kernel through nvcc, compile only (needs nvcc)
#   make bench         CUDA benchmark binary (needs nvcc)
#   make bench ARCH=sm_100
#   make device-roundtrip ARCH=sm_120
#                      folded tables built on a GPU, checked on the CPU
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
        pairtable.cuh fold_io.hpp

.PHONY: all test bench clean gen roundtrip tsan nvcc-check device-roundtrip
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

# The pair-table kernels themselves, compiled as host C++ behind a
# `__CUDACC__` shim and run with every emulated thread on a CPU thread of
# its own.  Everything above tests what the kernels are made of; this is
# the only thing that tests the kernels, and the fold against the CPU's
# own table.
test_emu_$(1): test_pairtable_emu.cpp $$(HDRS) curve_$(1).h vec_canon.h vec_fold.h
	$$(CXX) $$(CXXFLAGS) -pthread -DGPU_ECC2K_CURVE_HEADER='"curve_$(1).h"' -o $$@ $$<
endef
$(foreach c,$(CURVES),$(eval $(call TEST_RULE,$(c))))

test: $(TESTS)
	@rc=0; for t in $(TESTS); do \
	  echo "=== ./$$t ==="; ./$$t || rc=1; done; \
	echo "=== roundtrip ==="; $(MAKE) --no-print-directory roundtrip || rc=1; \
	$(MAKE) --no-print-directory tsan || rc=1; \
	exit $$rc

# The fold's round trip: the table the emulated kernel built, written out
# and taken over by `PairSumTable::from_folded_parts`, then asked every
# question the CPU's own table answers.  `test_emu_*` compares stored
# state; this compares answers, and is what runs every tag through the
# recovery that reads it.
#
# The plans come from `dump_fold_plan.rs` for the bases `vec_fold.h` was
# built from, and `test_emu_*` requires each to be `vec_fold.h`'s plan
# field by field -- which is what keeps the numbers below and the ones in
# `dump_fold_vectors.rs` from drifting apart.
FOLD_SEED    := 0x5eed_f01d
FOLD_BASE_k23 := 23 700
FOLD_BASE_k41 := 41 1300
FOLD_PLANS   := fold_plan_k23.bin fold_plan_k41.bin
FOLD_TABLES  := fold_table_k23.bin fold_table_k41.bin
fold_plan_%.bin: ../../examples/dump_fold_plan.rs ../../src/cryptanalysis/koblitz_fast.rs \
                 ../../src/cryptanalysis/koblitz_index_calculus.rs
	cd ../.. && $(CARGO) run --release --quiet --example dump_fold_plan -- \
	  $(FOLD_BASE_$*) $(FOLD_SEED) gpu/ecc2k/$@

fold_table_%.bin: test_emu_% fold_plan_%.bin
	./test_emu_$* --plan fold_plan_$*.bin --table $@ > /dev/null || { rm -f $@; exit 1; }

roundtrip: $(FOLD_TABLES)
	cd ../.. && $(CARGO) run --release --quiet --example load_fold_table -- \
	  $(addprefix gpu/ecc2k/,$(FOLD_TABLES))

# The same emulation under ThreadSanitizer.  The kernels share memory for
# real there -- the device build's counts, cursors and presence words --
# so a missing atomic or a shared scratch buffer is a report, not a
# timing accident.  And a control: `test_racy_k23` swaps each atomic for
# a plain read-modify-write, and has to be reported, in `pt_atomic_*`, or
# the clean runs would mean nothing.  It is reported in the count pass,
# and `halt_on_error` stops it there -- which matters beyond the exit
# code: a count that lost increments undersizes the word array, and a
# racy fill run past that point could write beyond it.
TSAN_FLAGS := -O1 -g -std=c++17 -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas \
              -fsanitize=thread -pthread
test_tsan_%: test_pairtable_emu.cpp $(HDRS) curve_%.h vec_canon.h vec_fold.h
	$(CXX) $(TSAN_FLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_$*.h"' -o $@ $<
test_racy_%: test_pairtable_emu.cpp $(HDRS) curve_%.h vec_canon.h vec_fold.h
	$(CXX) $(TSAN_FLAGS) -DPT_EMU_RACY -DGPU_ECC2K_CURVE_HEADER='"curve_$*.h"' -o $@ $<

tsan: test_tsan_k23 test_tsan_k41 test_racy_k23
	@rc=0; for c in k23 k41; do \
	  echo "=== ./test_tsan_$$c (ThreadSanitizer) ==="; \
	  if TSAN_OPTIONS=halt_on_error=1 ./test_tsan_$$c > tsan_$$c.log 2>&1; then \
	    tail -n 1 tsan_$$c.log; else cat tsan_$$c.log; rc=1; fi; \
	done; \
	echo "=== ./test_racy_k23 (every atomic a plain read-modify-write) ==="; \
	TSAN_OPTIONS=halt_on_error=1 ./test_racy_k23 > racy_k23.log 2>&1; \
	if grep -q "ThreadSanitizer: data race" racy_k23.log && grep -q "pt_atomic_" racy_k23.log; \
	then echo "reported as a data race in pt_atomic_*, as it must be"; \
	else cat racy_k23.log; echo "NOT reported: the ThreadSanitizer runs above prove nothing"; rc=1; fi; \
	exit $$rc

# Compile, don't run: every kernel here through nvcc.  Nothing else does --
# the host tests build the pair-table kernels with g++ behind a
# `__CUDACC__` shim, and `ptx_stats2k.sh` compiles the walk with clang --
# so this is the one check that the code a device would run gets past the
# compiler a device build would use.  No GPU needed.  Warnings are errors,
# and ptxas reports what each pair-table kernel takes in registers and
# stack.  `fold2k.cu` is built and linked, which covers `pairtable.cuh`
# and the launcher; `bench2k.cu` stands in for `kernels2k.cuh`.
NVCC_ARCHS ?= sm_90 sm_100 sm_120
NVCC_CHECK_FLAGS := -O3 -std=c++17 -Werror all-warnings
nvcc-check: $(HDRS) kernels2k.cuh bench2k.cu fold2k.cu curve_k23.h curve_k41.h curve_ecc2k95.h
	@for a in $(NVCC_ARCHS); do \
	  for c in k23 k41; do \
	    echo "=== nvcc -arch=$$a fold2k.cu, $$c ==="; \
	    out=$$($(NVCC) $(NVCC_CHECK_FLAGS) -arch=$$a -Xptxas -v -o nvcc_check_fold2k \
	      -DGPU_ECC2K_CURVE_HEADER="\"curve_$$c.h\"" fold2k.cu 2>&1) \
	      || { echo "$$out"; exit 1; }; \
	    echo "$$out" | grep -E "entry function|Used"; \
	  done; \
	  echo "=== nvcc -arch=$$a bench2k.cu, ecc2k95 ==="; \
	  $(NVCC) $(NVCC_CHECK_FLAGS) -arch=$$a -c -o /dev/null \
	    -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' bench2k.cu || exit 1; \
	done; rm -f nvcc_check_fold2k

# On a machine with a GPU: each fold table built on the device from the
# CPU's plan, then loaded and checked on the CPU.  The launcher is the
# only part of this the host tests cannot run.  For a wider base, write a
# plan with `dump_fold_plan.rs` directly and pass `--sample` to the
# loader (see README.md).
fold2k_%: fold2k.cu $(HDRS) curve_%.h
	$(NVCC) $(NVCCFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_$*.h"' -o $@ $<

DEVICE_TABLES := fold_table_dev_k23.bin fold_table_dev_k41.bin
device-roundtrip: fold2k_k23 fold2k_k41 $(FOLD_PLANS)
	./fold2k_k23 fold_plan_k23.bin fold_table_dev_k23.bin
	./fold2k_k41 fold_plan_k41.bin fold_table_dev_k41.bin
	cd ../.. && $(CARGO) run --release --quiet --example load_fold_table -- \
	  $(addprefix gpu/ecc2k/,$(DEVICE_TABLES))

bench: bench2k.cu $(HDRS) kernels2k.cuh curve_ecc2k95.h
	$(NVCC) $(NVCCFLAGS) -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' -o $@ $<

clean:
	rm -f $(GEN) $(TESTS) $(FOLD_PLANS) $(FOLD_TABLES) $(DEVICE_TABLES) bench \
	  test_tsan_k23 test_tsan_k41 test_racy_k23 tsan_*.log racy_*.log \
	  fold2k_k23 fold2k_k41 nvcc_check_fold2k
