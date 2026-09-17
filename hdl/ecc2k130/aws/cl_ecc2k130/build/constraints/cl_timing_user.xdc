# cl_timing_user.xdc -- timing constraints for cl_ecc2k130.
#
# Two clocks: clk_main_a0 (250 MHz, constrained by the shell's
# generated_cl_clocks_aws.xdc) carries the OCL port and the slow side of
# BRIDGE; clk_eng, which MMCM_ENG derives from it, carries everything
# else.  Vivado derives clk_eng's constraint from the MMCM attributes, so
# nothing here names a period; change the frequency in
# cl_ecc2k130_defines.vh (or CLK_MHZ in build_afi.sh), not here.
#
# The two clocks share a source, so Vivado would time paths between them
# against their edge relationship; the only such paths are inside BRIDGE
# (ec2k_axil_cdc.vhd), and they are the two kinds below.  Anything else
# crossing between the domains is a bug and is left timed so it fails
# loudly.
#
# Read at synthesis (PROCESSING_ORDER LATE), carried into implementation
# in the CL checkpoint.

# Flags: a level in one domain into an ASYNC_REG synchroniser in the other.
# The first stage is the only legitimate endpoint of a cross-domain path.
set_false_path -to [get_cells -hierarchical -filter {NAME =~ *BRIDGE*cdc_*_meta_reg*}]

# Data: written at least one source clock before its flag rises, sampled
# at least SYNC_FF destination clocks after, so it has a source period plus
# the synchroniser latency to settle; one source period is ample and keeps
# the routing honest.
set_max_delay -datapath_only 4.000 \
    -from [get_cells -hierarchical -filter {NAME =~ *BRIDGE*cdc_*_src_*_reg*}] \
    -to   [get_cells -hierarchical -filter {NAME =~ *BRIDGE*cdc_*_cap_*_reg*}]

# The engine-side reset comes from an xpm_cdc_async_rst, which carries its
# own constraints; MMCM_ENG/RST is not a timing endpoint.
