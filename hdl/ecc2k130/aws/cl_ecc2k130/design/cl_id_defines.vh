// cl_id_defines.vh -- PCIe identity of the ECC2K-130 image.
//
// Required by the AFI manifest.  Vendor 0x1D0F (Amazon) with a device id in
// the range 0xF000-0xF0FF the shell reserves for custom logic; the host
// program does not look at these, it recognises the image by the MAGIC
// register.

`ifndef CL_ID_DEFINES
`define CL_ID_DEFINES

  // 31:16 PCIe Device ID, 15:0 PCIe Vendor ID
  `define CL_SH_ID0       32'hF013_1D0F

  // 31:16 PCIe Subsystem ID, 15:0 PCIe Subsystem Vendor ID
  `define CL_SH_ID1       32'h2C13_FEDC

`endif
