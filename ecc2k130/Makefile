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
PACKED_POLY_CHAIN ?= 0
PACKED_UNROLL_INV ?= 0
PACKED_PAIR_PRODUCTS ?= 0
PACKED_POLY_STATE ?= 0
PACKED_DIRECT_REDUCE ?= 0
PACKED_GENERATED_PRODUCT ?= 0
PACKED_STATE_TILE ?= 0
# Keep experimental tuning independent: one change per comparison initially.
GPUFLAGS = -Xptxas -v
ifeq ($(GLOBAL_CG),1)
GPUFLAGS += -Xptxas --def-load-cache=cg
endif
# The CPU build picks the widest word the target supports (AVX-512 -> 512 lanes,
# AVX2 -> 256, otherwise 64).  Use a portable baseline when the build machine
# and the run machine may differ.
MARCH   ?= native
ARCH    ?= -gencode arch=compute_80,code=sm_80 \
           -gencode arch=compute_86,code=sm_86 \
           -gencode arch=compute_89,code=sm_89 \
           -gencode arch=compute_90,code=sm_90 \
           -gencode arch=compute_120,code=sm_120
# clang emits PTX for CLANGARCH; ptxas then compiles that PTX for SMARCH, which
# may be newer than any architecture clang itself knows about (e.g. sm_120).
CLANGARCH ?= sm_90
SMARCH  ?= sm_90

HDRS = include/bitslice.h include/fieldbs.h include/fieldpb.h include/walk.h \
       include/kernel.h include/ref.h include/solver.h include/bigmod.h \
       include/curveparams.h include/packed131.h include/packedtransform131.h include/packedsigma131.h include/packedpolyreduce131.h include/packeddirectreduce131.h include/packedgeneratedproduct131.h include/packedkernels.cuh include/packedengine.cuh
GEN  = generated/eccF131.h generated/eccF83.h generated/eccF41.h generated/eccF23.h \
       generated/eccP97.h generated/eccP41.h generated/eccP19.h generated/eccP13.h

DEFS = -DECC_BATCH=$(BATCH) -DECC_THREADS=$(THREADS) -DECC_MINBLOCKS=$(MINBLOCKS) \
       -DECC_STREAM_KARAT=$(STREAM_KARAT) -DECC_SMEM_SPILL=$(SMEM_SPILL) \
       -DECC_PACKED_SINGLE_PRODUCT=$(PACKED_SINGLE_PRODUCT) -DECC_PACKED_CACHE_DENOM=$(PACKED_CACHE_DENOM) \
       -DECC_PACKED_BY_VALUE=$(PACKED_BY_VALUE) -DECC_PACKED_PERM_SIGMA=$(PACKED_PERM_SIGMA) \
       -DECC_PACKED_POLY_CHAIN=$(PACKED_POLY_CHAIN) -DECC_PACKED_UNROLL_INV=$(PACKED_UNROLL_INV) \
       -DECC_PACKED_PAIR_PRODUCTS=$(PACKED_PAIR_PRODUCTS) -DECC_PACKED_POLY_STATE=$(PACKED_POLY_STATE) \
       -DECC_PACKED_DIRECT_REDUCE=$(PACKED_DIRECT_REDUCE) \
       -DECC_PACKED_GENERATED_PRODUCT=$(PACKED_GENERATED_PRODUCT) \
       -DECC_PACKED_STATE_TILE=$(PACKED_STATE_TILE)

all: cpu

# The generated leaf is sized to a register file: the largest halving-chain
# multiply whose peak liveness still fits REGS live values.  255 is a GPU's
# register budget and gives a 66-word leaf at m=131; a host has 16 general
# registers, so REGS=32 is the number to regenerate with when tuning the CPU
# reference build, and the two targets genuinely do not want the same leaf.
REGS ?= 255

generate:
	cd codegen && python3 gen.py --out ../generated --regs $(REGS)
	cd codegen && python3 genpacked.py
	cd codegen && python3 gensigma.py
	cd codegen && python3 genpolyreduce.py
	cd codegen && python3 gendirectreduce.py
	python3 codegen/genpackedproduct.py

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
	$(CXX) -O3 -std=c++17 -march=$(MARCH) -fopenmp $(DEFS) -DECC_NO_CUDA \
	    -Wall -Wno-unused-function -Wno-unknown-pragmas -x c++ src/main.cu -o $@

gpu: ecc2k130

ecc2k130: src/main.cu $(HDRS) $(GEN)
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) $(GPUFLAGS) -lineinfo -Xcompiler -O3 \
	    -Xcompiler -fopenmp src/main.cu -o $@ -lgomp

# Modal entry-point parameters must survive the trip through click; see
# codegen/checkcli.py.  Cheap, needs nothing installed, and the failure it
# catches is otherwise invisible until someone runs that entry point.
check-cli:
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

test-packed-cuda:
	mkdir -p build
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) src/testpackedcuda.cu -o build/test-packed-cuda
	./build/test-packed-cuda

test-timing:
	mkdir -p build
	$(CXX) -O1 -std=c++17 -Wno-unknown-pragmas src/testtiming.cpp -o build/test-timing
	./build/test-timing

test: check-cli check-cycles test-schedule test-packed test-packed-network test-generated-product test-timing check-verilog ecc2k130-cpu
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

# Offline device compile: no GPU needed, reports registers, spills and size.
ptx: src/ptxspike.cu $(HDRS) $(GEN)
	$(CLANG) -x cuda --cuda-device-only --cuda-gpu-arch=$(CLANGARCH) --cuda-path=$(CUDAPATH) \
	    -nocudalib -Wno-unknown-cuda-version -O3 -std=c++17 $(DEFS) -S -o build/spike.ptx src/ptxspike.cu
	$(PTXAS) -arch=$(SMARCH) -O3 -v build/spike.ptx -o build/spike.cubin

# Run on Modal GPUs (see modal_app.py).
modal:
	modal run modal_app.py::bench --gpu H100

# Audited packed configuration for the RTX PRO 6000 Blackwell Server Edition.
RTX_PRO6000_WORKERS ?= 192512
RTX_PRO6000_ENV = ECC_GPU=RTX-PRO-6000 ECC_CUDA_VERSION=13.3.1 \
                 ECC_PACKED_SINGLE_PRODUCT=1 ECC_PACKED_CACHE_DENOM=1 \
                 ECC_PACKED_BY_VALUE=1 ECC_PACKED_PERM_SIGMA=3 \
                 ECC_PACKED_POLY_CHAIN=1 ECC_PACKED_UNROLL_INV=1 \
                 ECC_PACKED_PAIR_PRODUCTS=1 ECC_PACKED_POLY_STATE=1 ECC_PACKED_DIRECT_REDUCE=1 \
                 ECC_PACKED_GENERATED_PRODUCT=1 ECC_PACKED_STATE_TILE=256

bench-rtx-pro6000:
	$(RTX_PRO6000_ENV) $(MODAL) run modal_app.py::bench --packed \
	    --batch 32 --threads 256 --leaf 0 --min-blocks 2 \
	    --steps 1024 --launches 32 --repeats 3 --workers $(RTX_PRO6000_WORKERS)

audit-rtx-pro6000:
	mkdir -p build
	$(RTX_PRO6000_ENV) $(MODAL) run packed_audit.py --min-blocks 2 --block-threads 256 --repeats 3 --workers $(RTX_PRO6000_WORKERS) \
	    --output build/rtx-pro6000-audit.json

clean:
	rm -f ecc2k130 ecc2k130-cpu build/*.ptx build/*.cubin build/*.o build/cudashim.h
	rm -f build/test-schedule-base build/test-schedule-stream
	rm -f build/test-packed build/test-packed-base build/test-packed-single
	rm -f build/test-packed-network build/test-packed-network-direct build/test-packed-network-generated
	rm -f build/test-packed-cuda
	rm -f build/test-timing

.PHONY: all generate check-cli check-cycles perfmodel test-decomp test-schedule test-packed test-packed-network test-generated-product test-packed-cuda test-timing rehearse autolab cpu gpu test break-small bench bench95 ptx check-cuda modal bench-rtx-pro6000 audit-rtx-pro6000 clean
