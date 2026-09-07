# hdl/ecc -- analyse, elaborate and run the self-checking testbenches.
#
#   make / make test   run every testbench
#   make vectors    regenerate vectors_secp256k1.txt from the Python oracle
#   make clean
#
# Requires GHDL (Debian/Ubuntu: apt install ghdl).  Everything is VHDL-2008
# with no vendor primitives, so it also reads cleanly into Vivado / Quartus.

GHDL  ?= ghdl
STD   ?= --std=08
PY    ?= python3

SRC := fp_pkg.vhd fp_mul_secp256k1.vhd ec_add_pipe.vhd
TBS := fp_mul_tb ec_add_tb
VEC := vectors_secp256k1.txt

.PHONY: all test vectors clean
all test: $(VEC)
	$(GHDL) -a $(STD) $(SRC) fp_mul_tb.vhd ec_add_tb.vhd
	@rc=0; for tb in $(TBS); do \
	  $(GHDL) -e $(STD) $$tb && \
	  echo "=== $$tb ===" && $(GHDL) -r $(STD) $$tb || rc=1; done; exit $$rc

vectors:
	$(PY) ../../gpu/ecc/ecref.py vhdl --curve secp256k1 > $(VEC)

$(VEC):
	$(PY) ../../gpu/ecc/ecref.py vhdl --curve secp256k1 > $(VEC)

clean:
	rm -f *.o *.cf $(TBS) e~*.o
