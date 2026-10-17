#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 9pt)
#text(size: 17pt, weight: "bold")[Held-out χ resultant solves: five-arm native subset]

*Run:* heldout-five-arm-v2 · *field:* `F_211^2` · *subgroup:* order 229 ·
*host:* Darwin arm64 shared L0.

#image("HELDOUT_GRAPH_V2.svg", width: 100%)

All 60 children returned verified scalars. A separate Rust checker
replayed both relation streams and all 332 relation witnesses.
χ and direct pair had identical trials, relations, rank-nine
recovery and scalars in all 12 paired rows. Mean cold S_group was
140.672 for χ, 88.908 for direct pair, 7.979 for signed rho and
4.361 for BSGS. The two-class floor is 0.886227.

The V1 setup ledger overcharged rho, BSGS and direct pair. V2
corrects it and retains V1. Source S3, quotient S3, direct
inverse-Vélu and full image on these targets are still pending;
portable primitive S and isolated wall are unset. Sources:
`HELDOUT_FIVE_ARM_SPEC_V2.json`, `HELDOUT_FIVE_ARM_V2.summary.json`,
`HELDOUT_FIVE_ARM_V2.check.json`, raw child files and hash manifest.
