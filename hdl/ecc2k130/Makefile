# hdl/ecc2k130 -- analyse, elaborate and run the self-checking testbenches.
#
#   make / make test   run every testbench
#   make check         prove the hardware algorithms against the client's field model
#   make vectors       regenerate vectors_ecc2k130.txt from that model
#   make cost          gate counts of the multiplier's constant matrices
#   make clean
#
# Requires GHDL (Debian/Ubuntu: apt install ghdl; macOS: the release
# tarball from github.com/ghdl/ghdl, the Homebrew cask is disabled).
# Everything is VHDL-2008 with no vendor primitives, so it also reads
# cleanly into Vivado / Quartus.

GHDL  ?= ghdl
STD   ?= --std=08
PY    ?= python3

SRC := gf131_pkg.vhd gf131_mul.vhd ec2k_step_pipe.vhd ec2k_walker.vhd
TBP := gf131_tb_pkg.vhd
TBS := gf131_mul_tb ec2k_step_tb ec2k_walker_tb
VEC := vectors_ecc2k130.txt

.PHONY: all test check vectors cost clean
all test: $(VEC)
	$(GHDL) -a $(STD) $(SRC) $(TBP) $(addsuffix .vhd,$(TBS))
	@rc=0; for tb in $(TBS); do \
	  $(GHDL) -e $(STD) $$tb && \
	  echo "=== $$tb ===" && $(GHDL) -r $(STD) $$tb || rc=1; done; exit $$rc

check:
	$(PY) ecc2k_ref.py check

vectors:
	$(PY) ecc2k_ref.py vectors > $(VEC)

$(VEC):
	$(PY) ecc2k_ref.py vectors > $(VEC)

cost:
	$(PY) ecc2k_ref.py cost

clean:
	rm -f *.o *.cf $(TBS) e~*.o
