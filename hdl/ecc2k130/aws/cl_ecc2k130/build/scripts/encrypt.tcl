# encrypt.tcl -- gather the sources into build/src_post_encryption.
#
# The HDK's stock encrypt.tcl copies only Verilog from ${design_dir}.  The
# engine is VHDL two directories up (hdl/ecc2k130/*.vhd), so this version
# copies the RTL list from there as well, then encrypts both languages the
# way the stock script does.  build_all.tcl sources this file from
# build/scripts with design_dir, src_post_enc_dir, HDK_SHELL_DIR,
# HDK_SHELL_DESIGN_DIR and ENCRYPT already set.

if {[llength [glob -nocomplain -dir ${src_post_enc_dir} *]] != 0} {
  eval file delete -force [glob ${src_post_enc_dir}/*]
}

# shell interface includes
foreach f [glob -directory ${HDK_SHELL_DESIGN_DIR}/interfaces *.inc] {
  file copy -force $f ${src_post_enc_dir}/
}

# the CL top and its defines
foreach f [glob -directory ${design_dir} *.{v,sv,vh,svh,inc}] {
  file copy -force $f ${src_post_enc_dir}/
}

# the engine: RTL only, testbenches stay behind.  ECC_RTL_DIR overrides the
# location for a CL directory that has been copied out of the repository.
if {[info exists ::env(ECC_RTL_DIR)]} {
  set rtl_dir $::env(ECC_RTL_DIR)
} else {
  set rtl_dir [file normalize ${design_dir}/../../..]
}
set rtl_files {
  gf131_pkg.vhd
  gf2_kmul.vhd
  gf131_mul.vhd
  ec2k_batch_pipe.vhd
  ec2k_walker.vhd
  ec2k_axil.vhd
}
foreach f $rtl_files {
  if {![file exists ${rtl_dir}/$f]} {
    puts "ERROR: engine source ${rtl_dir}/$f not found"
    exit 1
  }
  file copy -force ${rtl_dir}/$f ${src_post_enc_dir}/
}

exec chmod +w {*}[glob ${src_post_enc_dir}/*]

if {$ENCRYPT} {
  print "Encryption enabled. Encrypting HDL files and DCPs."
  encrypt -k ${HDK_SHELL_DIR}/build/scripts/vivado_keyfile.txt      -lang verilog -quiet [glob -nocomplain -- ${src_post_enc_dir}/*.{v,sv,vh,inc}]
  encrypt -k ${HDK_SHELL_DIR}/build/scripts/vivado_vhdl_keyfile.txt -lang vhdl    -quiet [glob -nocomplain -- ${src_post_enc_dir}/*.vhd?]
} else {
  print "Encryption disabled."
}
