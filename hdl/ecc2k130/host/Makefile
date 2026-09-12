# hdl/ecc2k130/host -- the FPGA host program.
#
#   make               ecc2k130-fpga with the software model only (any machine)
#   make pci           ecc2k130-fpga with the F2 PCI backend as well
#                      (needs the aws-fpga SDK: source sdk_setup.sh first, or
#                      set SDK_DIR)
#   make test          selftest, then a --sim run at a loose dp weight through
#                      the whole worker.py contract: dp file, checkpoint,
#                      resume, exit codes
#
# The program is one translation unit over the client's headers in
# ../../../ecc2k130/include, built exactly like the client's own
# test-packed target.

CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wno-unknown-pragmas -Wno-unused-function
SDK_DIR  ?= $(shell echo $$SDK_DIR)

BIN := ecc2k130-fpga
SRC := ecc2k130_fpga.cpp
INC := ../../../ecc2k130/include

.PHONY: all pci test clean
all: $(BIN)

$(BIN): $(SRC) $(wildcard $(INC)/*.h)
	$(CXX) $(CXXFLAGS) $(SRC) -o $@

pci: $(SRC) $(wildcard $(INC)/*.h)
	@test -n "$(SDK_DIR)" || { echo "SDK_DIR is not set: source aws-fpga/sdk_setup.sh"; exit 1; }
	$(CXX) $(CXXFLAGS) -DECC_FPGA_PCI -I$(SDK_DIR)/userspace/include $(SRC) \
	    -L$(SDK_DIR)/userspace/lib/so -lfpga_mgmt -Wl,-rpath,$(SDK_DIR)/userspace/lib/so -o $(BIN)

test: $(BIN)
	./test.sh ./$(BIN)

clean:
	rm -f $(BIN)
