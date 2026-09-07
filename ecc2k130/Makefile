# ECC2K-130 client.
#
#   make generate     regenerate the arithmetic headers from codegen/
#   make cpu          CPU build (64-bit lanes, OpenMP)
#   make test         run the validation suite
#   make gpu          CUDA build (needs nvcc)
#   make ptx          device-only compile + ptxas report (needs clang or nvcc,
#                     but no GPU) -- register, spill and code-size numbers
#
# Build-time knobs:
#   BATCH      walks batched into one field inversion per thread (default 32)
#   THREADS    CUDA block size (default 128)
#   MINBLOCKS  resident blocks per SM to ask ptxas for (default 2)
#   ARCH     CUDA architectures

NVCC    ?= nvcc
CXX     ?= g++
CLANG   ?= clang
PTXAS   ?= ptxas
CUDAPATH ?= /usr/local/cuda
BATCH   ?= 32
THREADS ?= 128
MINBLOCKS ?= 2
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
       include/curveparams.h
GEN  = generated/eccF131.h generated/eccF83.h generated/eccF41.h generated/eccF23.h \
       generated/eccP97.h generated/eccP41.h generated/eccP19.h generated/eccP13.h

DEFS = -DECC_BATCH=$(BATCH) -DECC_THREADS=$(THREADS) -DECC_MINBLOCKS=$(MINBLOCKS)

all: cpu

# The generated leaf is sized to a register file: the largest halving-chain
# multiply whose peak liveness still fits REGS live values.  255 is a GPU's
# register budget and gives a 66-word leaf at m=131; a host has 16 general
# registers, so REGS=32 is the number to regenerate with when tuning the CPU
# reference build, and the two targets genuinely do not want the same leaf.
REGS ?= 255

generate:
	cd codegen && python3 gen.py --out ../generated --regs $(REGS)

cpu: ecc2k130-cpu

ecc2k130-cpu: src/main.cu $(HDRS) $(GEN)
	$(CXX) -O3 -std=c++17 -march=$(MARCH) -fopenmp $(DEFS) -DECC_NO_CUDA \
	    -Wall -Wno-unused-function -Wno-unknown-pragmas -x c++ src/main.cu -o $@

gpu: ecc2k130

ecc2k130: src/main.cu $(HDRS) $(GEN)
	$(NVCC) -O3 -std=c++17 $(ARCH) $(DEFS) -lineinfo -Xcompiler -O3 \
	    -Xcompiler -fopenmp src/main.cu -o $@ -lgomp

# Modal entry-point parameters must survive the trip through click; see
# codegen/checkcli.py.  Cheap, needs nothing installed, and the failure it
# catches is otherwise invisible until someone runs that entry point.
check-cli:
	python3 codegen/checkcli.py modal_app.py

test: check-cli ecc2k130-cpu
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

clean:
	rm -f ecc2k130 ecc2k130-cpu build/*.ptx build/*.cubin build/*.o build/cudashim.h

.PHONY: all generate check-cli cpu gpu test break-small bench bench95 ptx check-cuda modal clean
