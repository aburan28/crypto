// cl_id_defines.vh -- PCIe identity of the ECC2K-130 image.
//
// Required by the AFI manifest.  Vendor 0x1D0F (Amazon) with a device id in
// the range 0xF000-0xF0FF the shell reserves for custom logic; the host
// program does not look at these, it recognises the image by the MAGIC
// register.
//
// The HDK's aws_build_dcp_from_cl.py copies these into the manifest with
// str.lstrip("32'h"), so an id whose first hex digit is 2 or 3 loses it
// (the first image had subsystem id 2C13, the manifest said C13, and the
// loaded FPGA was refused with cl-id-mismatch).  build_afi_instance.sh now
// checks the manifest against this file; the ids here also avoid the trap.

`ifndef CL_ID_DEFINES
`define CL_ID_DEFINES

  // 31:16 PCIe Device ID, 15:0 PCIe Vendor ID
  `define CL_SH_ID0       32'hF013_1D0F

  // 31:16 PCIe Subsystem ID, 15:0 PCIe Subsystem Vendor ID
  `define CL_SH_ID1       32'hEC13_FEDC

`endif
