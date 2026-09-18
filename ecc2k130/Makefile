# ECC2K-130 client.
#
#   make generate     regenerate the arithmetic headers from codegen/
#   make cpu          CPU build (64-bit lanes, OpenMP)
#   make test         run the validation suite
#   make gpu          CUDA build (needs nvcc)
#   make ptx          device-only compile + ptxas report (needs clang or nvcc,
#                     but no GPU) -- register, spill and code-size numbers
#   make verilog      emit the FPGA datapath from the same IR (python3 only)
#   make check-verilog  simulate the emitted RTL against the field model
#
# Build-time knobs:
#   BATCH      walks batched into one field inversion per thread (default 32)
#   THREADS    CUDA block size (default 128)
#   MINBLOCKS  resident blocks per SM to ask ptxas for (default 2)
#   PACKED_SINGLE_PRODUCT  use one polynomial product for packed mul (default 0)
#   ARCH     CUDA architectures

MODAL   ?= modal
NVCC    ?= nvcc
CXX     ?= g++
# Apple's clang needs the runtime named separately; everywhere else the driver
# flag is enough.  Only the host walk is threaded, so this is all it reaches.
#   make cpu CXX=clang++ OMPFLAGS="-Xpreprocessor -fopenmp -lomp"
OMPFLAGS ?= -fopenmp
CLANG   ?= clang
PTXAS   ?= ptxas
PTX     ?= /tmp/autolab.ptx
CUDAPATH ?= /usr/local/cuda
BATCH   ?= 32
THREADS ?= 128
MINBLOCKS ?= 2
STREAM_KARAT ?= 0
SMEM_SPILL ?= 0
GLOBAL_CG ?= 0
PACKED_SINGLE_PRODUCT ?= 0
PACKED_CACHE_DENOM ?= 0
PACKED_BY_VALUE ?= 0
PACKED_PERM_SIGMA ?= 0
PACKED_PARTIAL_SIGMA ?= 0
PACKED_SIGMA_ORDER ?= 0
PACKED_BYTE_SIGMA ?= 0
PACKED_FAST_CONVERT ?= 0
PACKED_POLY_CHAIN ?= 0
PACKED_UNROLL_INV ?= 0
PACKED_PAIR_PRODUCTS ?= 0
PACKED_POLY_STATE ?= 0
PACKED_DIRECT_REDUCE ?= 0
PACKED_GENERATED_PRODUCT ?= 0
PACKED_CLMAD ?= 0
PACKED_NATIVE_PRODUCT ?= 0
PACKED_NATIVE_REDUCE ?= 0
PACKED_FUSED_SIGMA ?= 0
PACKED_INLINE ?= 0
PACKED_WEIGHTED_PREFIX ?= 0
PACKED_COMPACT_STATE ?= 0
PACKED_SHARED_SIGMA ?= 0
PACKED_STATE_TILE ?= 0
PACKED_BLOCK_INVERSE ?= 0
PACKED_BATCH_SPLIT ?= 1
PACKED_LAST_SLOT_CACHE ?= 0
PACKED_TAIL_LAYOUT ?= 0
PACKED_SHARED_X_SLOTS ?= 0
PACKED_XONLY_23 ?= 0
PACKED_XONLY_DOUBLE_ONLY ?= 0
PACKED_XONLY_BRIDGE3 ?= 0
PACKED_XONLY_BRIDGE1_COMMON ?= 0
PACKED_XONLY_BRIDGE_MOD72 ?= 0
PACKED_XONLY_SKIP_EMPTY_BRIDGE ?= 0
PACKED_XONLY_ARITHMETIC_ONLY ?= 0
PACKED_XONLY_POLY_SELECT ?= 0
PACKED_XONLY_POLY_DP_CONVERT ?= 1
PACKED_XONLY_POLY_SELECT_MASK ?= 4095
PACKED_XONLY_POLY_SELECT_TARGET ?= 14
# Keep experimental tuning independent: one change per comparison initially.
GPUFLAGS = -Xptxas -v
ifeq ($(GLOBAL_CG),1)
GPUFLAGS += -Xptxas --def-load-cache=cg
endif
# The CPU build picks the widest word the target supports (AVX-512 -> 512 lanes,
# AVX2 -> 256, NEON -> 128, otherwise 64).  Use a portable baseline when the
# build machine and the run machine may differ.  That lever does not work on
# Apple clang, whose driver keeps the host CPU's features regardless of -march:
# MARCH=armv8-a leaves __ARM_FEATURE_SHA3 undefined yet the backend still folds
# XOR trees into EOR3.  Build portable binaries on the target OS.
MARCH   ?= native

# On aarch64 the wide word is NEON, and EOR3 (FEAT_SHA3) is most of what makes
# it faster than the 64-bit fallback rather than merely wider.  -march=native
# does not always imply +sha3 -- GCC on Darwin reports no SHA3 for native on
# cores that have it -- so ask the compiler which spelling actually defines the
# feature macro and take the first that does.  A compiler-accepted +sha3
# spelling is not proof the host can execute EOR3 (armv8.2-a+sha3 is accepted
# on any recent toolchain), so those candidates are used only when the CPU
# reports FEAT_SHA3.  They also rewrite the requested -march, so only search
# when MARCH is native: an explicit MARCH names a target this host is not
# evidence about.  Falling through to plain MARCH still gets NEON and BSL,
# just not the three-input XOR.
HOSTARCH := $(shell uname -m)
HOSTOS   := $(shell uname -s)
ifneq (,$(filter arm64 aarch64,$(HOSTARCH)))
# Linux lists sha3 in cpuinfo Features when FEAT_SHA3 is present.  Darwin has
# no /proc, and every Apple Silicon core has SHA3, so a missing
# hw.optional.arm.FEAT_SHA3 key means present there and absent anywhere else.
ifeq ($(HOSTOS),Darwin)
HOSTSHA3 := $(shell sysctl -n hw.optional.arm.FEAT_SHA3 2>/dev/null || echo 1)
else
HOSTSHA3 := $(shell grep -m1 '^Features' /proc/cpuinfo 2>/dev/null | grep -c -w sha3)
endif
sha3flag = $(shell $(CXX) $(1) -dM -E -x c++ /dev/null 2>/dev/null | grep -c __ARM_FEATURE_SHA3)
SHA3TRY := $(if $(filter native,$(MARCH)), \
    $(if $(filter-out 0,$(HOSTSHA3)), \
        -march=native+sha3 -mcpu=native+sha3 -march=armv8.2-a+sha3))
ARCHFLAGS := $(firstword \
    $(foreach f,-march=$(MARCH) $(SHA3TRY), \
        $(if $(filter-out 0,$(call sha3flag,$(f))),$(f))) \
    -march=$(MARCH))
else
ARCHFLAGS := -march=$(MARCH)
endif
ARCH    ?= -gencode arch=compute_80,code=sm_80 \
           -gencode arch=compute_86,code=sm_86 \
           -gencode arch=compute_89,code=sm_89 \
           -gencode arch=compute_90,code=sm_90 \
           -gencode arch=compute_120,code=sm_120
# clang emits PTX for CLANGARCH; ptxas then compiles that PTX for SMARCH, which
# may be newer than any architecture clang itself knows about (e.g. sm_120).
CLANGARCH ?= sm_90
SMARCH  ?= sm_90

HDRS = include/durable.h include/bitslice.h include/fieldbs.h include/fieldpb.h include/walk.h \
       include/kernel.h include/ref.h include/solver.h include/bigmod.h \
       include/curveparams.h include/packed131.h include/packednative131.h include/packedfusedsigma131.h include/packedtransform131.h include/packedsigma131.h include/packedpolyreduce131.h include/packeddirectreduce131.h include/packedgeneratedproduct131.h include/packedcompactstate.cuh include/packedblockinverse131.cuh include/packedblockinversereference131.cuh include/packedxonly23.cuh include/packedkernels.cuh include/packedengine.cuh
GEN  = generated/eccF131.h generated/eccF83.h generated/eccF41.h generated/eccF23.h \
       generated/eccP97.h generated/eccP41.h generated/eccP19.h generated/eccP13.h

DEFS = -DECC_BATCH=$(BATCH) -DECC_THREADS=$(THREADS) -DECC_MINBLOCKS=$(MINBLOCKS) \
       -DECC_STREAM_KARAT=$(STREAM_KARAT) -DECC_SMEM_SPILL=$(SMEM_SPILL) \
       -DECC_PACKED_SINGLE_PRODUCT=$(PACKED_SINGLE_PRODUCT) -DECC_PACKED_CACHE_DENOM=$(PACKED_CACHE_DENOM) \
       -DECC_PACKED_BY_VALUE=$(PACKED_BY_VALUE) -DECC_PACKED_PERM_SIGMA=$(PACKED_PERM_SIGMA) \
       -DECC_PACKED_PARTIAL_SIGMA=$(PACKED_PARTIAL_SIGMA) \
       -DECC_PACKED_SIGMA_ORDER=$(PACKED_SIGMA_ORDER) \
       -DECC_PACKED_BYTE_SIGMA=$(PACKED_BYTE_SIGMA) \
       -DECC_PACKED_FAST_CONVERT=$(PACKED_FAST_CONVERT) \
       -DECC_PACKED_POLY_CHAIN=$(PACKED_POLY_CHAIN) -DECC_PACKED_UNROLL_INV=$(PACKED_UNROLL_INV) \
       -DECC_PACKED_PAIR_PRODUCTS=$(PACKED_PAIR_PRODUCTS) -DECC_PACKED_POLY_STATE=$(PACKED_POLY_STATE) \
       -DECC_PACKED_DIRECT_REDUCE=$(PACKED_DIRECT_REDUCE) \
       -DECC_PACKED_GENERATED_PRODUCT=$(PACKED_GENERATED_PRODUCT) \
       -DECC_PACKED_CLMAD=$(PACKED_CLMAD) \
       -DECC_PACKED_NATIVE_PRODUCT=$(PACKED_NATIVE_PRODUCT) \
       -DECC_PACKED_NATIVE_REDUCE=$(PACKED_NATIVE_REDUCE) \
       -DECC_PACKED_FUSED_SIGMA=$(PACKED_FUSED_SIGMA) \
       -DECC_PACKED_INLINE=$(PACKED_INLINE) \
       -DECC_PACKED_WEIGHTED_PREFIX=$(PACKED_WEIGHTED_PREFIX) \
       -DECC_PACKED_COMPACT_STATE=$(PACKED_COMPACT_STATE) \
       -DECC_PACKED_SHARED_SIGMA=$(PACKED_SHARED_SIGMA) \
       -DECC_PACKED_STATE_TILE=$(PACKED_STATE_TILE) \
       -DECC_PACKED_BLOCK_INVERSE=$(PACKED_BLOCK_INVERSE) \
       -DECC_PACKED_BATCH_SPLIT=$(PACKED_BATCH_SPLIT) \
       -DECC_PACKED_LAST_SLOT_CACHE=$(PACKED_LAST_SLOT_CACHE) \
       -DECC_PACKED_TAIL_LAYOUT=$(PACKED_TAIL_LAYOUT) \
       -DECC_PACKED_SHARED_X_SLOTS=$(PACKED_SHARED_X_SLOTS) \
       -DECC_PACKED_XONLY_23=$(PACKED_XONLY_23) \
       -DECC_PACKED_XONLY_DOUBLE_ONLY=$(PACKED_XONLY_DOUBLE_ONLY) \
       -DECC_PACKED_XONLY_BRIDGE3=$(PACKED_XONLY_BRIDGE3) \
       -DECC_PACKED_XONLY_BRIDGE1_COMMON=$(PACKED_XONLY_BRIDGE1_COMMON) \
       -DECC_PACKED_XONLY_BRIDGE_MOD72=$(PACKED_XONLY_BRIDGE_MOD72) \
       -DECC_PACKED_XONLY_SKIP_EMPTY_BRIDGE=$(PACKED_XONLY_SKIP_EMPTY_BRIDGE) \
       -DECC_PACKED_XONLY_ARITHMETIC_ONLY=$(PACKED_XONLY_ARITHMETIC_ONLY) \
       -DECC_PACKED_XONLY_POLY_SELECT=$(PACKED_XONLY_POLY_SELECT) \
       -DECC_PACKED_XONLY_POLY_DP_CONVERT=$(PACKED_XONLY_POLY_DP_CONVERT) \
       -DECC_PACKED_XONLY_POLY_SELECT_MASK=$(PACKED_XONLY_POLY_SELECT_MASK) \
       -DECC_PACKED_XONLY_POLY_SELECT_TARGET=$(PACKED_XONLY_POLY_SELECT_TARGET)

all: cpu

# The generated leaf is sized to a register file: the largest halving-chain
# multiply whose peak liveness still fits REGS live values.  255 is a GPU's
# register budget and gives a 66-word leaf at m=131.
#
# A smaller leaf for the host is the obvious thing to want and it is measurably
# wrong.  Single core, M4 Pro, NEON words, curve 131, median of five:
#
#   REGS=32   leaf 9     9.820 M it/s
#   REGS=64   leaf 17   10.052
#   REGS=128  leaf 33   10.799
#   REGS=255  leaf 66   11.353
#
# The leaf spills either way -- 59% of the emitted multiply is loads and stores
# even at leaf 66 -- and spilling to a hot stack is cheaper than the extra
# Karatsuba level a small leaf buys.  Both targets want the same leaf.
REGS ?= 255

generate:
	cd codegen && python3 gen.py --out ../generated --regs $(REGS)
	cd codegen && python3 genpacked.py
	cd codegen && python3 gensigma.py
	cd codegen && python3 genpolyreduce.py
	cd codegen && python3 gendirectreduce.py
	python3 codegen/genpackedproduct.py
	python3 codegen/genfusedsigma.py

# FPGA datapath: the same IR emitted as pipelined Verilog instead of C.
# See ../FPGA-CEILING.md; needs python3 only, no FPGA toolchain.
verilog:
	cd codegen && python3 genverilog.py --out ../generated/rtl

# Read the emitted RTL back and simulate it against field.Onb.  Checks the
# emission (stage crossings, register chains, output wiring), not the IR.
check-verilog:
	cd codegen && python3 testverilog.py --rtl ../generated/rtl

cpu: ecc2k130-cpu

ecc2k130-cpu: src/main.cu $(HDRS) $(GEN)
	$(CXX) -O3 -std=c++17 $(ARCHFLAGS) $(OMPFLAGS) $(DEFS) -DECC_NO_CUDA \
	    -Wall -Wno-unused-function -Wno-unknown-pragmas -x c++ src/main.cu -o $@

gpu: ecc2k130

ecc2k130: src/main.cu $(HDRS) $(GEN)
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) $(GPUFLAGS) -lineinfo -Xcompiler -O3 \
	    -Xcompiler -fopenmp src/main.cu -o $@ -lgomp

# Modal entry-point parameters must survive the trip through click; see
# codegen/checkcli.py.  Cheap, needs nothing installed, and the failure it
# catches is otherwise invisible until someone runs that entry point.
check-cli:
	python3 -m unittest discover -s codegen -p testbenchmarkhardware.py
	python3 codegen/genfusedsigma.py --check
	python3 codegen/genpacked.py --check --out include/packedtransform131.h
	python3 codegen/gensigma.py --check --out include/packedsigma131.h
	python3 codegen/genpolyreduce.py --check --out include/packedpolyreduce131.h
	python3 codegen/gendirectreduce.py --check --out include/packeddirectreduce131.h
	python3 codegen/checkcli.py modal_app.py
	python3 codegen/checkcli.py profile_diagnostic.py
	python3 codegen/checkcli.py packed_audit.py
	python3 codegen/checkcli.py polystate_audit.py
	python3 codegen/autolab.py --self-test
	python3 -m unittest discover -s codegen -p 'test*report.py'
	python3 -m unittest discover -s codegen -p 'testautotune*.py'
	python3 -m unittest discover -s codegen -p 'testpolystateaudit.py'
	python3 -m unittest discover -s codegen -p 'testdirectbuild.py'

# The iteration exploits an automorphism group of order 2m, and the usual price
# for that is fruitless cycles.  This iteration cannot have them -- it is
# equivariant, so it needs no canonical representative and has nothing for a
# sign to feed back through -- but the arithmetic coincidence that would still
# produce one is checkable, so check it, for every curve, on every change to the
# step schedule.
check-cycles:
	cd codegen && python3 cycles.py --max-length 8

# SAT-based point decomposition (codegen/decomp.py, indexcalc.py).  Kept out of
# `test` because it needs solvers the rest of the repo does not:
#   pip install python-sat pycryptosat
test-decomp:
	cd codegen && python3 testdecomp.py

# Drive the whole campaign shape on a curve that finishes in seconds: several
# workers with disjoint seed spaces, one resumed from its checkpoint, and a
# collision that only the merge can see because no single worker holds both
# halves.  Those are the parts a long run only reaches at the very end.
rehearse: ecc2k130-cpu
	./rehearse.sh

build/test-production: src/testproduction.cpp src/main.cu $(HDRS) $(GEN)
	mkdir -p build
	$(CXX) -O2 -std=c++17 -DECC_NO_CUDA -DECC_BATCH=2 -Wno-unknown-pragmas -Wno-unused-function src/testproduction.cpp -o $@

test-production: ecc2k130-cpu build/test-production
	./build/test-production
	python3 -m unittest discover -s aws -p test_certification.py -v

certify-local: ecc2k130-cpu build/test-production
	python3 aws/certify.py --output build/certification-local.json

# Offline search over the build space against ptxas output, which needs no GPU.
# It ranks nothing it can measure the speed of -- it produces the few
# configurations worth spending GPU time on, for `modal run
# modal_app.py::autotune` to actually measure.  Uses nvcc when CUDAPATH has one
# and the check-cuda clang route when it does not; `modal run
# modal_app.py::autolab` runs the same search in a CPU container that has nvcc.
# Where a walk iteration actually goes, weighted by how many times a step calls
# each routine, plus the model's error against every rate ever measured.  Needs
# PTX for the build in question; make ptx or autolab leaves one behind.
perfmodel:
	cd codegen && python3 perfmodel.py --ptx $(PTX)

autolab:
	cd codegen && python3 autolab.py --cuda-path=$(CUDAPATH) --ptxas=$(PTXAS) \
	    --clang=$(CLANG) --arch=$(patsubst sm_%,%,$(SMARCH))

test-schedule:
	mkdir -p build
	$(CXX) -O2 -std=c++17 -DECC_STREAM_KARAT=0 -Wno-unknown-pragmas src/testschedule.cpp -o build/test-schedule-base
	$(CXX) -O2 -std=c++17 -DECC_STREAM_KARAT=1 -Wno-unknown-pragmas src/testschedule.cpp -o build/test-schedule-stream
	./build/test-schedule-base
	./build/test-schedule-stream

test-packed:
	mkdir -p build
	$(CXX) -O2 -std=c++17 -DECC_PACKED_SINGLE_PRODUCT=0 -DECC_PACKED_DIRECT_REDUCE=0 -DECC_PACKED_GENERATED_PRODUCT=0 -Wno-unknown-pragmas src/testpacked.cpp -o build/test-packed-base
	$(CXX) -O2 -std=c++17 -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_DIRECT_REDUCE=0 -DECC_PACKED_GENERATED_PRODUCT=0 -Wno-unknown-pragmas src/testpacked.cpp -o build/test-packed-single
	./build/test-packed-base
	./build/test-packed-single

test-packed-network:
	mkdir -p build
	$(CXX) -O2 -std=c++17 -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_BY_VALUE=1 -DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_DIRECT_REDUCE=0 -DECC_PACKED_GENERATED_PRODUCT=0 -Wno-unknown-pragmas src/testpacked.cpp -o build/test-packed-network
	$(CXX) -O2 -std=c++17 -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_BY_VALUE=1 -DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_GENERATED_PRODUCT=0 -Wno-unknown-pragmas src/testpacked.cpp -o build/test-packed-network-direct
	$(CXX) -O2 -std=c++17 -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_BY_VALUE=1 -DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_GENERATED_PRODUCT=1 -Wno-unknown-pragmas src/testpacked.cpp -o build/test-packed-network-generated
	./build/test-packed-network
	./build/test-packed-network-direct
	./build/test-packed-network-generated

test-generated-product:
	python3 codegen/genpackedproduct.py --check
	python3 -m unittest discover -s codegen -p testgeneratedproduct.py

test-clmad:
	python3 -m unittest discover -s codegen -p testclmad.py
	mkdir -p build
	$(CXX) -O2 -std=c++17 -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_BY_VALUE=1 -DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_GENERATED_PRODUCT=1 -DECC_PACKED_CLMAD=1 -Wno-unknown-pragmas src/testpacked.cpp -o build/test-packed-clmad-host
	./build/test-packed-clmad-host

test-packed-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) src/testpackedcuda.cu -o build/test-packed-cuda
	./build/test-packed-cuda

test-packed-block-inverse-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) src/testblockinversecuda.cu -o build/test-packed-block-inverse-cuda
	./build/test-packed-block-inverse-cuda

test-packed-storage-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) src/testpackedstatecuda.cu -o build/test-packed-storage-cuda
	./build/test-packed-storage-cuda

test-shared-sigma-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) src/testsharedsigmacuda.cu -o build/test-shared-sigma-cuda
	./build/test-shared-sigma-cuda

test-timing:
	mkdir -p build
	$(CXX) -O1 -std=c++17 -Wno-unknown-pragmas src/testtiming.cpp -o build/test-timing
	./build/test-timing

# Every host word type this machine can build, against the scalar definitions.
# Built at the tuned architecture so the wide word under test is the one the
# walk will actually use.
test-logic:
	mkdir -p build
	$(CXX) -O2 -std=c++17 $(ARCHFLAGS) -Wno-unknown-pragmas src/testlogic.cpp -o build/test-logic
	./build/test-logic

test: check-cli check-cycles test-schedule test-packed test-packed-network test-generated-product test-clmad test-timing test-logic check-verilog ecc2k130-cpu
	./ecc2k130-cpu --test

# End-to-end: recover planted discrete logarithms on the small curves, through
# both field backends.
break-small: ecc2k130-cpu
	@echo "-- normal basis --"
	@for i in 0 1 2 3; do \
	    ./ecc2k130-cpu --curve 23 --instance $$i --threads 1 --steps 32 --verify 2 | grep -E "planted|MISMATCH"; \
	done
	@for i in 0 1 2 3; do \
	    ./ecc2k130-cpu --curve 41 --instance $$i --threads 4 --steps 64 --verify 2 | grep -E "planted|MISMATCH"; \
	done
	@echo "-- polynomial basis --"
	@for i in 0 1 2 3; do \
	    ./ecc2k130-cpu --curve 19 --instance $$i --threads 1 --steps 16 --verify 2 | grep -E "planted|MISMATCH"; \
	done
	@for i in 0 1 2 3; do \
	    ./ecc2k130-cpu --curve 41 --poly-basis --instance $$i --threads 4 --steps 32 --verify 2 | grep -E "planted|MISMATCH"; \
	done

# Throughput on ECC2K-95's curve.
bench95: ecc2k130-cpu
	./ecc2k130-cpu --curve 97 --bench --steps 32 --launches 8 --verify 0

bench: ecc2k130-cpu
	./ecc2k130-cpu --curve 131 --bench --steps 32 --launches 8 --verify 0

# Type-check the whole CUDA translation unit, host and device, with clang.
# Needs no GPU and no nvcc, only a CUDA include tree at CUDAPATH.  On a machine
# with no toolkit installed, pip supplies one:
#   pip install nvidia-cuda-runtime-cu12 nvidia-cuda-nvcc-cu12 \
#               nvidia-curand-cu12 nvidia-cuda-cccl-cu12
# then point CUDAPATH at a directory holding those include trees merged into
# include/, and PTXAS at nvidia/cuda_nvcc/bin/ptxas.  CLANGARCH must be an
# architecture the local clang knows: clang 18 stops at sm_90.
check-cuda: src/main.cu $(HDRS) $(GEN)
	printf '#include <cuda_runtime.h>\nextern "C" cudaError_t cudaConfigureCall(dim3, dim3, size_t, cudaStream_t);\n' > build/cudashim.h
	$(CLANG) -x cuda --cuda-host-only --cuda-path=$(CUDAPATH) -nocudalib \
	    -Wno-unknown-cuda-version -O2 -std=c++17 $(DEFS) -include build/cudashim.h \
	    -c src/main.cu -o build/main_host.o
	$(CLANG) -x cuda --cuda-device-only --cuda-gpu-arch=$(CLANGARCH) --cuda-path=$(CUDAPATH) \
	    -nocudalib -Wno-unknown-cuda-version -O3 -std=c++17 $(DEFS) \
	    -include build/cudashim.h -S -o build/main_dev.ptx src/main.cu
	$(PTXAS) -arch=$(SMARCH) -O3 -v build/main_dev.ptx -o build/main_dev.cubin

# Where the packed walk's instructions go: no GPU needed.  Prints the ALU
# composition of a scalar update and what a throughput objective implies.
# See THROUGHPUT-30B.md.
path-cost:
	CUDAPATH=$(CUDAPATH) CLANG=$(CLANG)++ python3 path_cost.py $(PATH_COST_ARGS)

# Offline device compile: no GPU needed, reports registers, spills and size.
ptx: src/ptxspike.cu $(HDRS) $(GEN)
	$(CLANG) -x cuda --cuda-device-only --cuda-gpu-arch=$(CLANGARCH) --cuda-path=$(CUDAPATH) \
	    -nocudalib -Wno-unknown-cuda-version -O3 -std=c++17 $(DEFS) -S -o build/spike.ptx src/ptxspike.cu
	$(PTXAS) -arch=$(SMARCH) -O3 -v build/spike.ptx -o build/spike.cubin

# Run on Modal GPUs (see modal_app.py).
modal:
	modal run modal_app.py::bench --gpu H100

# Audited packed configuration for the RTX PRO 6000 Blackwell Server Edition.
RTX_PRO6000_WORKERS ?= 385024
RTX_PRO6000_COMPACT_STATE ?= 1
RTX_PRO6000_SHARED_SIGMA ?= 1
RTX_PRO6000_ENV = ECC_GPU=RTX-PRO-6000 ECC_CUDA_VERSION=13.3.1 \
                 ECC_PACKED_SINGLE_PRODUCT=1 ECC_PACKED_CACHE_DENOM=1 \
                 ECC_PACKED_BY_VALUE=1 ECC_PACKED_PERM_SIGMA=3 \
                 ECC_PACKED_POLY_CHAIN=1 ECC_PACKED_UNROLL_INV=1 \
                 ECC_PACKED_PAIR_PRODUCTS=1 ECC_PACKED_POLY_STATE=1 ECC_PACKED_DIRECT_REDUCE=1 \
                 ECC_PACKED_GENERATED_PRODUCT=1 ECC_PACKED_CLMAD=1 ECC_PACKED_STATE_TILE=256 \
                 ECC_PACKED_WEIGHTED_PREFIX=2 ECC_PACKED_COMPACT_STATE=$(RTX_PRO6000_COMPACT_STATE) \
                 ECC_PACKED_SHARED_SIGMA=$(RTX_PRO6000_SHARED_SIGMA)

bench-rtx-pro6000:
	$(RTX_PRO6000_ENV) $(MODAL) run modal_app.py::bench --packed \
	    --batch 16 --threads 256 --leaf 0 --min-blocks 2 \
	    --steps 1024 --launches 32 --repeats 3 --workers $(RTX_PRO6000_WORKERS)

audit-rtx-pro6000:
	mkdir -p build
	$(RTX_PRO6000_ENV) $(MODAL) run packed_audit.py --batch 16 --min-blocks 2 --block-threads 256 --repeats 3 --workers $(RTX_PRO6000_WORKERS) \
	    --output build/rtx-pro6000-audit.json

clean:
	rm -f ecc2k130 ecc2k130-cpu build/*.ptx build/*.cubin build/*.o build/cudashim.h
	rm -f build/test-schedule-base build/test-schedule-stream
	rm -f build/test-packed build/test-packed-base build/test-packed-single
	rm -f build/test-packed-network build/test-packed-network-direct build/test-packed-network-generated
	rm -f build/test-packed-clmad-host
	rm -f build/test-packed-cuda build/test-packed-storage-cuda build/test-packed-block-inverse-cuda
	rm -f build/test-shared-sigma-cuda build/test-last-slot-cache-cuda build/test-tail-product-cuda build/test-shared-x-cache-cuda
	rm -f build/test-timing build/test-logic

.PHONY: path-cost all generate check-cli check-cycles perfmodel test-decomp test-schedule test-packed test-packed-network test-generated-product test-clmad test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda test-timing test-logic test-production certify-local rehearse autolab cpu gpu test break-small bench bench95 ptx check-cuda modal bench-rtx-pro6000 audit-rtx-pro6000 clean

# Local packed preset. Prefer the measured optional assembler-only prefix,
# then pinned CUDA 13.3, then NVCC. Native CLMAD needs CUDA 13.3 or newer.
# Separate output and forced rebuild prevent stale binaries when knobs change.
LOCAL_GPU_ARCH ?= -gencode arch=compute_120,code=sm_120
LOCAL_PACKED_NVCC ?= $(if $(wildcard build/assembler134-toolchain/bin/nvcc),build/assembler134-toolchain/bin/nvcc,$(if $(wildcard build/cuda133/nvidia/cu13/bin/nvcc),build/cuda133/nvidia/cu13/bin/nvcc,$(NVCC)))
LOCAL_PACKED_CLMAD ?= $(shell $(LOCAL_PACKED_NVCC) --version 2>/dev/null | awk '/release / { split($$5,v,"."); print (v[1]+0>13 || (v[1]+0==13 && v[2]+0>=3)); found=1 } END { if (!found) print 0 }')
LOCAL_PACKED_INLINE ?= $(if $(filter 1,$(LOCAL_PACKED_CLMAD)),3,0)
# Exact seven-stage basis conversion, confirmed on the native G7 preset.
LOCAL_PACKED_FAST_CONVERT ?= $(if $(filter 1,$(LOCAL_PACKED_CLMAD)),1,0)
# Native normal products avoid basis conversions inside inversion. The
# software fallback keeps its established single-product implementation.
LOCAL_PACKED_SINGLE_PRODUCT ?= $(if $(filter 1,$(LOCAL_PACKED_CLMAD)),0,1)
# Share inversion across eight full warps in the validated native preset.
LOCAL_PACKED_BLOCK_INVERSE ?= $(if $(filter 1,$(LOCAL_PACKED_CLMAD)),1,0)
LOCAL_PACKED_BATCH ?= 16
# Exact 131-port walk routing is measured on the native RTX PRO 4500 preset.
LOCAL_PACKED_PARTIAL_SIGMA ?= $(if $(filter 1,$(LOCAL_PACKED_CLMAD)),1,0)
# Measured Frobenius instruction reductions in the native G7 preset.
LOCAL_PACKED_SIGMA_ORDER ?= $(if $(and $(filter 1,$(LOCAL_PACKED_CLMAD)),$(filter 1 3,$(LOCAL_PACKED_PARTIAL_SIGMA))),1,0)
LOCAL_PACKED_BYTE_SIGMA ?= $(if $(and $(filter 1,$(LOCAL_PACKED_CLMAD)),$(filter 1 3,$(LOCAL_PACKED_PARTIAL_SIGMA))),1,0)
# The measured RTX PRO 4500 native preset uses two physical workers per batch.
LOCAL_PACKED_BATCH_SPLIT ?= $(if $(and $(filter 1,$(LOCAL_PACKED_CLMAD)),$(filter 1,$(LOCAL_PACKED_BLOCK_INVERSE)),$(filter 16,$(LOCAL_PACKED_BATCH))),2,1)
LOCAL_PACKED_MINBLOCKS ?= $(if $(filter 2,$(LOCAL_PACKED_BATCH_SPLIT)),3,2)
LOCAL_PACKED_LAST_SLOT_CACHE ?= $(if $(and $(filter 1,$(LOCAL_PACKED_CLMAD)),$(filter 1,$(LOCAL_PACKED_BLOCK_INVERSE)),$(filter 16,$(LOCAL_PACKED_BATCH)),$(filter 2,$(LOCAL_PACKED_BATCH_SPLIT))),1,0)
# Regroup the high three product bits using measured native funnel shifts.
LOCAL_PACKED_TAIL_LAYOUT ?= $(if $(filter 1,$(LOCAL_PACKED_CLMAD)),1,0)
LOCAL_PACKED_SHARED_X_SLOTS ?= $(if $(and $(filter 1,$(LOCAL_PACKED_CLMAD)),$(filter 1,$(LOCAL_PACKED_BLOCK_INVERSE)),$(filter 16,$(LOCAL_PACKED_BATCH)),$(filter 2,$(LOCAL_PACKED_BATCH_SPLIT))),4,0)
LOCAL_PACKED_FLAGS = NVCC="$(LOCAL_PACKED_NVCC)" ARCH="$(LOCAL_GPU_ARCH)" \
    BATCH=$(LOCAL_PACKED_BATCH) THREADS=256 MINBLOCKS=$(LOCAL_PACKED_MINBLOCKS) \
    PACKED_SINGLE_PRODUCT=$(LOCAL_PACKED_SINGLE_PRODUCT) PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 \
    PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
    PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 \
    PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=$(LOCAL_PACKED_CLMAD) PACKED_INLINE=$(LOCAL_PACKED_INLINE) PACKED_WEIGHTED_PREFIX=2 \
    PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 PACKED_STATE_TILE=256 \
    PACKED_BLOCK_INVERSE=$(LOCAL_PACKED_BLOCK_INVERSE) PACKED_BATCH_SPLIT=$(LOCAL_PACKED_BATCH_SPLIT) PACKED_PARTIAL_SIGMA=$(LOCAL_PACKED_PARTIAL_SIGMA) \
    PACKED_SIGMA_ORDER=$(LOCAL_PACKED_SIGMA_ORDER) PACKED_BYTE_SIGMA=$(LOCAL_PACKED_BYTE_SIGMA) \
    PACKED_FAST_CONVERT=$(LOCAL_PACKED_FAST_CONVERT) \
    PACKED_LAST_SLOT_CACHE=$(LOCAL_PACKED_LAST_SLOT_CACHE) PACKED_TAIL_LAYOUT=$(LOCAL_PACKED_TAIL_LAYOUT) \
    PACKED_SHARED_X_SLOTS=$(LOCAL_PACKED_SHARED_X_SLOTS) \
    PACKED_XONLY_23=1 PACKED_XONLY_DOUBLE_ONLY=1 PACKED_XONLY_BRIDGE3=1 \
    PACKED_XONLY_BRIDGE1_COMMON=1 PACKED_XONLY_BRIDGE_MOD72=1 \
    PACKED_XONLY_SKIP_EMPTY_BRIDGE=$(PACKED_XONLY_SKIP_EMPTY_BRIDGE)
LOCAL_PACKED_WORKERS ?= 524288
LOCAL_PACKED_STEPS ?= 1024
LOCAL_PACKED_LAUNCHES ?= 16

# Explicit opt-in install; the regular build never downloads a compiler.
setup-local-cuda133:
	python3 -m pip install --no-deps --upgrade --target build/cuda133 \
	    nvidia-cuda-nvcc==13.3.73 nvidia-cuda-runtime==13.3.29 \
	    nvidia-nvvm==13.3.73 nvidia-cuda-crt==13.3.73

# Explicit opt-in: retain the 13.3 front end and replace only ptxas.
setup-local-assembler134: setup-local-cuda133
	python3 setup_cuda_assembler.py

.PHONY: setup-local-cuda133 setup-local-assembler134 gpu-local-packed bench-local-packed

gpu-local-packed:
	$(MAKE) -B build/ecc2k130-local-packed $(LOCAL_PACKED_FLAGS)

build/ecc2k130-local-packed: src/main.cu $(HDRS) $(GEN)
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) $(GPUFLAGS) -lineinfo \
	    -Xcompiler -O3 -Xcompiler -fopenmp src/main.cu -o $@ -lgomp

bench-local-packed: gpu-local-packed
	python3 benchmark_hardware.py
	./build/ecc2k130-local-packed --packed --bench --curve 131 \
	    --threads $(LOCAL_PACKED_WORKERS) --steps $(LOCAL_PACKED_STEPS) \
	    --launches $(LOCAL_PACKED_LAUNCHES) --verify 0

test-last-slot-cache-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) $(GPUFLAGS) -lineinfo \
	    -Xcompiler -O3 -Xcompiler -fopenmp src/testlastslotcuda.cu -o build/test-last-slot-cache-cuda -lgomp
	./build/test-last-slot-cache-cuda

test-tail-product-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) $(GPUFLAGS) -lineinfo \
	    -Xcompiler -O3 -Xcompiler -fopenmp src/testtailproductcuda.cu -o build/test-tail-product-cuda -lgomp
	./build/test-tail-product-cuda

test-shared-x-cache-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) $(GPUFLAGS) -lineinfo \
	    -Xcompiler -O3 -Xcompiler -fopenmp src/testsharedxcuda.cu -o build/test-shared-x-cache-cuda -lgomp
	./build/test-shared-x-cache-cuda

test-local-packed:
	$(MAKE) test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda $(if $(filter 1,$(LOCAL_PACKED_BLOCK_INVERSE)),test-packed-block-inverse-cuda) $(if $(filter 1 2,$(LOCAL_PACKED_LAST_SLOT_CACHE)),test-last-slot-cache-cuda) $(if $(filter 1,$(LOCAL_PACKED_CLMAD)),test-tail-product-cuda) $(if $(filter 2 4,$(LOCAL_PACKED_SHARED_X_SLOTS)),test-shared-x-cache-cuda) $(LOCAL_PACKED_FLAGS)

.PHONY: setup-local-cuda133 gpu-local-packed bench-local-packed test-local-packed test-packed-block-inverse-cuda test-last-slot-cache-cuda test-tail-product-cuda

.PHONY: test-shared-x-cache-cuda
