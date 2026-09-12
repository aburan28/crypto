# synth_cl_ecc2k130.tcl -- synthesis of the ECC2K-130 custom logic, in the
# form the F2 HDK's build_all.tcl expects (it sources this after
# encrypt.tcl, then build_level_1_cl.tcl for implementation).
#
# Differences from synth_CL_TEMPLATE.tcl: the engine is VHDL-2008, read
# with read_vhdl; and the geometry defines come from the environment
# (ECC_NENG, ECC_ID_W, ECC_DP_WEIGHT, ...), which build_afi.sh sets, so one
# source tree builds every image.

# Common header
source ${HDK_SHELL_DIR}/build/scripts/synth_cl_header.tcl


###############################################################################
print "Reading user source codes"
###############################################################################

# the CL top (SystemVerilog) ...
read_verilog -sv [glob ${src_post_enc_dir}/*.{s,}v]

# ... and the engine (VHDL-2008), all in the default library
read_vhdl -vhdl2008 [glob ${src_post_enc_dir}/*.vhd]

###############################################################################
print "Reading user constraints"
###############################################################################

read_xdc [ list \
  ${constraints_dir}/cl_synth_user.xdc \
  ${constraints_dir}/cl_timing_user.xdc
]

set_property PROCESSING_ORDER LATE [get_files cl_synth_user.xdc]
set_property PROCESSING_ORDER LATE [get_files cl_timing_user.xdc]

###############################################################################
print "Geometry"
###############################################################################

# Every ECC_* define present in the environment becomes a -verilog_define;
# cl_ecc2k130_defines.vh supplies the defaults for the rest.
set geometry {}
foreach name {ECC_NENG ECC_ID_W ECC_LOG_W ECC_LOG_NB ECC_FLUSH_CLK ECC_DP_WEIGHT ECC_DP_FIFO_W} {
  if {[info exists ::env($name)]} {
    lappend geometry -verilog_define ${name}=$::env($name)
    print "  $name = $::env($name)"
  }
}
if {[llength $geometry] == 0} {
  print "  (defaults from cl_ecc2k130_defines.vh)"
}

###############################################################################
print "Starting synthesizing customer design ${CL}"
###############################################################################
update_compile_order -fileset sources_1

synth_design -mode out_of_context \
             -top ${CL} \
             -verilog_define XSDB_SLV_DIS \
             {*}$geometry \
             -part ${DEVICE_TYPE} \
             -keep_equivalent_registers

# The numbers the capacity estimate in ../../../README.md waits for.
report_utilization -file ${reports_dir}/${CL}.${TAG}.synth_utilization.rpt
# Depth 4 reaches inside one engine: walker, step unit, multiplier.
report_utilization -hierarchical -hierarchical_depth 4 \
    -file ${reports_dir}/${CL}.${TAG}.synth_utilization_hier.rpt
report_timing_summary -max_paths 20 \
    -file ${reports_dir}/${CL}.${TAG}.synth_timing.rpt


# Common footer
source ${HDK_SHELL_DIR}/build/scripts/synth_cl_footer.tcl
