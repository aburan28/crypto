# synth.tcl -- out-of-context synthesis of the secp256k1 datapath, to turn
# the estimates in docs/ecc_fpga_cost_model.md into measurements.
#
#   vivado -mode batch -source synth.tcl -tclargs <part> <period_ns> <outdir>
#
# Out-of-context means no I/O buffers and no board constraints: the design
# is synthesised as a black box, which is what you want for a datapath that
# will later be wrapped in the AWS shell.  The numbers it produces are
# post-synthesis, not post-route, so treat the timing as an upper bound --
# placement and routing on a part as large as a VU47P typically costs
# another 10-20% of the clock.  Utilisation, by contrast, is close to final
# for a DSP-bound design like this one.
#
# Each module is synthesised separately so the DSP and LUT cost of the
# multiplier can be read off on its own, which is the number the cost model
# actually depends on.

if {$argc < 3} {
    puts "usage: vivado -mode batch -source synth.tcl -tclargs <part> <period_ns> <outdir>"
    exit 1
}
set part    [lindex $argv 0]
set period  [lindex $argv 1]
set outdir  [lindex $argv 2]
set rtldir  [file normalize [file join [file dirname [info script]] ..]]

file mkdir $outdir
puts "== part      $part"
puts "== period    $period ns ([format %.1f [expr {1000.0/$period}]] MHz)"
puts "== rtl from  $rtldir"

# Each entry: top-level entity, and the sources it needs.
set targets {
    fp_mul_secp256k1 {fp_pkg.vhd fp_mul_secp256k1.vhd}
    ec_add_pipe      {fp_pkg.vhd fp_mul_secp256k1.vhd ec_add_pipe.vhd}
}

set summary {}

foreach {top srcs} $targets {
    puts "\n======== synthesising $top ========"
    create_project -in_memory -part $part
    foreach s $srcs {
        read_vhdl -vhdl2008 [file join $rtldir $s]
    }
    # -mode out_of_context keeps the tool from inserting I/O buffers on
    # what are internal ports in the real design.
    synth_design -top $top -part $part -mode out_of_context

    # Constrain after synthesis: the clock is the only constraint this
    # datapath has, and every path in it is register-to-register.
    create_clock -period $period -name clk [get_ports clk]
    set_input_delay  -clock clk 0.0 [remove_from_collection [all_inputs] [get_ports clk]]
    set_output_delay -clock clk 0.0 [all_outputs]

    report_utilization -file [file join $outdir ${top}_util.rpt]
    report_utilization -hierarchical -file [file join $outdir ${top}_util_hier.rpt]
    report_timing_summary -max_paths 10 \
        -file [file join $outdir ${top}_timing.rpt]
    write_checkpoint -force [file join $outdir ${top}_synth.dcp]

    # Pull the headline numbers back out for a one-line summary.  Counting
    # cells by primitive group is more robust across Vivado versions than
    # scraping report_utilization's text.
    set wns [get_property SLACK [lindex [get_timing_paths -delay_type max -max_paths 1] 0]]
    if {$wns eq "" || $wns eq "NONE"} { set wns 0 }
    set fmax [expr {1000.0 / ($period - $wns)}]

    set nlut [llength [get_cells -hier -filter {PRIMITIVE_GROUP == LUT}]]
    set nff  [llength [get_cells -hier -filter {PRIMITIVE_GROUP == FLOP_LATCH}]]
    set ndsp [llength [get_cells -hier -filter {PRIMITIVE_GROUP == ARITHMETIC}]]
    set nram [llength [get_cells -hier -filter {PRIMITIVE_GROUP == BLOCKRAM}]]

    lappend summary [format "%-20s LUT %7d  FF %7d  DSP %5d  BRAM %4d  WNS %+7.3f ns  Fmax %6.1f MHz" \
        $top $nlut $nff $ndsp $nram $wns $fmax]
    close_project
}

set fh [open [file join $outdir summary.txt] w]
puts "\n================ SUMMARY ================"
puts $fh "part $part, target period $period ns"
foreach line $summary {
    puts $line
    puts $fh $line
}
puts $fh ""
puts $fh "Fmax is post-synthesis: place-and-route typically costs another 10-20%."
close $fh
puts "\nreports written to $outdir"
