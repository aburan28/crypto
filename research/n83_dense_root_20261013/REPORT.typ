#set page(width: 297mm, height: 210mm, margin: (x: 12mm, y: 9mm))
#set text(font: "Helvetica Neue", size: 8pt)
#set heading(numbering: none)

= N83 packed-root capacity and prefix-base result
K1798 dense · 298,468 usable points · 1,798 folded columns · 8 GiB allocated root slots

K1798 dense is the lowest-mean-probe completed prefix base in the frozen 64-point comparison: 169,069.8125 probes per target, versus 300,021.90625 for K1500 dense. The paired mean difference is 130,952.09375 probes, 95% bootstrap interval [52,057.703125, 209,990.828125]. Complete archived-base, four-point relation and scalar replay passed for all 256 target rows in the four completed arms.

#figure(image("reports/run-20261013T0942Z/primary.svg", width: 100%), caption: [Single public target: exact probes, root-table allocation, and exploratory online time. Source run: run-20261013T0942Z.])

The same public target was solved in 63.324125 ms online by K1798 dense and 1,173.29325 ms online by the selected strong-rho control, giving an exploratory rho/IC ratio of 18.528377 on the shared Apple M4 Pro host. Both recovered scalar 2927171373054142 and passed independent generic-curve replay. IC preparation was 319.165 s setup, 19.825 s rank and 0.014 s final matrix algebra, outside its target-online timer. K1500 regular reached exit 124 at the 900 s cap in both primary and panel runs.

The exact curve is degree-83 Koblitz $a=b=1$ over polynomial $z^83+z^7+z^4+z^2+1$, subgroup order 8,569,786,107,849,059. All IC arms use rank seed 29, 12 rank threads, one target thread and the same producer binary. K1271 regular/dense first witnesses matched on all 64 points; dense mode reduced its allocated packed slots from 8 GiB to 4 GiB. Six Lean 4 current-capacity theorems built without axioms; two conditional 75%-load theorems also built, one using standard `propext`. The Rust-loop bridge remains partial. The proposed 75% rule is not a measured arm.

Evidence: `analysis_primary.json`, `analysis_panel.json`, `runs/`, `PROTOCOL.md`, `PROOF_MAP.md`, and the two run reports. CPU isolation was not established on the shared host; controlled speedup remains null.

#pagebreak()

= Frozen panel and proof

#figure(image("reports/run-20261013T1018Z/panel.svg", width: 90%), caption: [Mean target probes, allocated slots and exploratory online time on the identical 64 public points. Source run: run-20261013T1018Z.])

#figure(image("reports/run-20261013T1018Z/panel_paired.svg", width: 70%), caption: [Paired mean-probe differences, left minus right, with 10,000-resample 95% bootstrap intervals. Positive favors the right-hand arm.])
