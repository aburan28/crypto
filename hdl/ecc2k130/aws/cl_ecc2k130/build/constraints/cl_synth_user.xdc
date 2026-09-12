# cl_synth_user.xdc -- synthesis-time constraints for cl_ecc2k130.
#
# None needed: the design is one clock domain (clk_main_a0, constrained by
# the shell's generated_cl_clocks_aws.xdc) and uses no vendor primitives.
# The file exists because synth_cl_ecc2k130.tcl reads it, as the template
# does.
