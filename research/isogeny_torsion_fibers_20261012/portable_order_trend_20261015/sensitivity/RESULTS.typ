#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)

= How large is the missing rho price?

*Frozen native derivation, 2026-10-15.* Four complete-solve count files,
their independent replay certificates and 48 paired prefix/round rows
passed hash, key, scalar-count and group-total checks. No solver input
or child was changed. Host: Darwin arm64 shared L0; elapsed time is not
used as a performance result.

For each child, let $G_("full")$ and $G_("rho")$ count cold group adds
and doubles from validation through final verification, and let
$R_("steps")$ count all observed rho walk steps including retries.

$W_("break") = (G_("full") - G_("rho")) / R_("steps")$.

This is the *extra omitted rho charge per observed step* needed to tie
the full-image group-count lower bound if every omitted full-image
operation is priced at zero. A positive omitted full-image price makes
the required rho charge larger. It is not a measured native conversion,
complete portable $S$, or end-to-end speed ratio.

#table(
  columns: (0.6fr, 1.5fr, 1fr, 1fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 5pt,
  [*Order $r$*], [*N64 median GAE/step*], [*Round min*], [*Round max*],
  [229], [112.973], [99.105], [131.619],
  [421], [303.573], [302.425], [351.781],
  [433], [334.844], [297.947], [344.670],
  [457], [288.595], [267.223], [299.550],
)

#image("SENSITIVITY_V1.svg", width: 100%)

Across every prefix and round, the minimum is 99.105 GAE/step at N64,
order 229; the maximum is 46,594.667 at N1, order 421. The N1
threshold charges the cold full-image census to one solve.

The source count data, exact run hashes and replay certificates are in
`../TREND_V1.json`, `PROTOCOL_V1.md` and `RUN.md`; all 48 calculated
rows are in `SENSITIVITY_V1.json`. The native source is
`examples/different_j_price_sensitivity.rs`. The parent report retains
the matched BSGS and rho rows, separate Fp multiplication/inversion
axes and memory. Complete portable S still requires disjoint native
operation counts, calibration across host classes and independent
trace coverage; isolated wall and second-host performance are separate
gates. A direct chi-coordinate oracle that avoids the source census is
the structural next test.
