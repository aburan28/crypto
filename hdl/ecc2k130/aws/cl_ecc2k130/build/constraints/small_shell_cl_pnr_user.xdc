# small_shell_cl_pnr_user.xdc -- user floorplan for cl_ecc2k130 with the F2
# small shell.  build_level_1_cl.tcl reads ${SHELL_MODE}_cl_pnr_user.xdc
# after the shell's own floorplan (small_shell_level_1_fp_cl.xdc, which
# defines pblock_CL, the whole dynamic region across the three SLRs), so
# the file has to exist even when it says nothing.
#
# It says nothing on purpose.  The design is NENG identical, independent
# walker engines (WRAPPER/CL/ENGINE/engines[*].eng) hanging off one small
# AXI-Lite register block, with no DDR, HBM or PCIM traffic; the only
# inter-engine nets are the load/report buses to ENGINE.  The placer keeps
# each engine's ~15k LUTs together on its own and chooses the SLR per
# engine; a hand floorplan would only pin how many engines go to each SLR,
# which the HDK examples do because they have DDR controllers and DMA
# crossbars to keep apart.  If a build shows SLR-crossing paths inside an
# engine in the post-route timing report, add child pblocks of pblock_CL
# here (see hdk/cl/examples/*/build/constraints/small_shell_cl_pnr_user.xdc
# for the clock-region ranges of each SLR) and split engines[*] across them
# roughly 2:1:1 for SLR2:SLR1:SLR0, which is the usable area ratio.
