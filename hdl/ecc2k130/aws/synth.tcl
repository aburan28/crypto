# synth.tcl -- out-of-context synthesis of the ECC2K-130 engine, to replace
# the LUT and clock estimates in ../README.md with Vivado's numbers.
#
#   vivado -mode batch -source synth.tcl -tclargs <part> <period_ns> <outdir>
#
# Out-of-context means no I/O buffers and no board constraints: the design
# is synthesised as a black box, which is what you want for a datapath that
# will later sit inside the AWS shell.  Numbers are post-synthesis; treat
# the clock as an upper bound (place-and-route on a VU47P typically costs
# another 10-20%).  Utilisation is close to final: this design is all LUTs
# and flip-flops -- binary-field arithmetic uses no DSP48 at all -- and LUT
# counts move little between synthesis and routing.
#
# Each module is synthesised separately so the cost of one multiplier, the
# batched step unit, the simple step unit and the whole walker can be read
# off on their own.  The multiplier's LUT count is the number the per-part
# capacity estimate in ../README.md depends on; the multiplier is also the
# place to sweep MUL_KARATSUBA (gf131_pkg.vhd) for the LUT optimum.
#
# ../../ecc/aws/run_aws_synthesis.sh launches a build instance and runs this
# file; point it here with HDLDIR=$(pwd)/.. -- see README.md alongside.

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
set mul {gf131_pkg.vhd gf2_kmul.vhd gf131_mul.vhd}
set targets [list \
    gf131_mul       $mul \
    ec2k_batch_pipe [concat $mul ec2k_batch_pipe.vhd] \
    ec2k_step_pipe  [concat $mul ec2k_step_pipe.vhd] \
    ec2k_walker     [concat $mul ec2k_batch_pipe.vhd ec2k_walker.vhd] \
]

set summary {}

foreach {top srcs} $targets {
    puts "\n======== synthesising $top ========"
    create_project -in_memory -part $part
    foreach s $srcs {
        read_vhdl -vhdl2008 [file join $rtldir $s]
    }
    synth_design -top $top -part $part -mode out_of_context

    # The clock is the only constraint; every path is register-to-register.
    create_clock -period $period -name clk [get_ports clk]
    set_input_delay  -clock clk 0.0 [remove_from_collection [all_inputs] [get_ports clk]]
    set_output_delay -clock clk 0.0 [all_outputs]

    report_utilization -file [file join $outdir ${top}_util.rpt]
    report_utilization -hierarchical -file [file join $outdir ${top}_util_hier.rpt]
    report_timing_summary -max_paths 10 \
        -file [file join $outdir ${top}_timing.rpt]
    write_checkpoint -force [file join $outdir ${top}_synth.dcp]

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
