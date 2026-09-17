# cl_synth_user.xdc -- synthesis-time constraints for cl_ecc2k130.
#
# None needed.  The two clocks are constrained by the shell's
# generated_cl_clocks_aws.xdc (clk_main_a0) and derived from MMCM_ENG
# (clk_eng); the crossings between them are in cl_timing_user.xdc.  The
# file exists because synth_cl_ecc2k130.tcl reads it, as the template does.
