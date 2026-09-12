# cl_timing_user.xdc -- timing constraints for cl_ecc2k130.
#
# Every path in the engine is register to register on clk_main_a0, which
# the shell constrains at 250 MHz.  The reset pipeline is the only thing to
# say anything about: rst_pipe is a two-flop synchroniser of a signal that
# is already synchronous, so it needs no exception, and ec2k_axil holds
# every engine in reset through a registered CTRL bit rather than the
# shell reset directly.
#
# If the first build misses timing, do not relax it here; see
# ../../README.md for the two places in the datapath to look.
