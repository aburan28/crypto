# hdl/ecc2k130 -- the rho sequencer around the generated ECC2K-130 datapath.
#
#   make / make test   build and run every self-checking testbench
#   make clean
#
# Requires Icarus Verilog (brew install icarus-verilog, apt install iverilog).
# The datapath modules and the reference vectors are generated, so run
#
#   make -C ../../ecc2k130 verilog
#
# first.  Everything here is plain Verilog-2001 with no vendor primitives, so it
# also reads into Vivado or Quartus unchanged.
#
# The cost model these belong to is ../../ecc2k130/FPGA-CEILING.md.

IVERILOG ?= iverilog
VVP      ?= vvp
STD      ?= -g2012
GEN      ?= ../../ecc2k130/generated/rtl

GENSRC := $(GEN)/ecc_mul131.v $(GEN)/ecc_delay.v $(GEN)/ecc_inv131.v \
          $(GEN)/ecc_pre131.v $(GEN)/ecc_post131.v
SRC    := ecc_rho_core.v
TBS    := tb_ecc_rho_core

.PHONY: all test clean
all test:
	@rc=0; for tb in $(TBS); do \
	  $(IVERILOG) $(STD) -o $$tb.vvp $(GENSRC) $(SRC) $$tb.v || exit 1; \
	  echo "=== $$tb ==="; $(VVP) $$tb.vvp || rc=1; done; exit $$rc

clean:
	rm -f *.vvp
