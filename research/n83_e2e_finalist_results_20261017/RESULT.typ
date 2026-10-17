#set page(width: 297mm, height: 210mm, margin: (x: 15mm, y: 13mm))
#set text(font: "Helvetica Neue", size: 9pt)
#set heading(numbering: none)

= N83 factor-base frontier on sixteen fresh Tier 1 points

K1107/root11 lowered one-target cold process cost by *3.949 seconds* against
K1120/root11 on the same sixteen public points (paired 95% interval
[3.404, 4.494] seconds, exact sign-flip $p = 0.000030518$). K1107/root22
lowered it by 4.332 seconds [3.318, 5.346]. Both K1107 layouts beat every
other tested K in cold paired mean under the frozen Bonferroni gate. Their
mutual comparison remained unresolved, so the eight-configuration cold
selection is `BOUNDED_UNKNOWN`; the online selection is also
`BOUNDED_UNKNOWN`.

#table(
  columns: (33mm, 26mm, 25mm, 29mm, 27mm, 33mm, 28mm),
  inset: 4pt,
  stroke: 0.4pt + luma(200),
  [*Configuration*], [*Usable points*], [*Columns*], [*Online median ms*],
  [*Cold mean s*], [*Median rho / IC*], [*IC faster / 16*],
  [K1104/root11], [183,264], [1,104], [113.394], [108.382], [16.247], [16],
  [K1104/root22], [183,264], [1,104], [111.651], [108.112], [11.334], [16],
  [K1107/root11], [183,762], [1,107], [79.820], [104.826], [18.507], [16],
  [K1107/root22], [183,762], [1,107], [82.929], [104.635], [23.910], [16],
  [K1112/root11], [184,592], [1,112], [195.899], [106.596], [12.457], [15],
  [K1112/root22], [184,592], [1,112], [191.117], [106.699], [15.636], [16],
  [K1120/root11], [185,920], [1,120], [73.487], [108.775], [22.355], [16],
  [K1120/root22], [185,920], [1,120], [71.913], [108.966], [19.394], [16],
)

Every arm ran in a fresh process with an empty base, index and log cache.
The online interval starts with target-dependent work after that reusable
preparation and ends with scalar verification. The cold process interval
charges the preparation. Each arm has an independently verified rho solve
of the same point. Two distinct GitHub-hosted AMD EPYC 9V74 allocations,
each exposing four vCPUs, passed Tier 1 topology, resource and before/after
rho calibration gates. There are 128 IC–rho pairs, four IC repeat controls
and 256 rho calibrations. The point, base, four-point relation, scalar,
source, binary and resource records passed independent local replay.

The source commit is `416489e09c3d8885b6c12f262935e6ba29373460`;
the IC source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.
The raw, provider-bound receipts are retained in a 23,533,732-byte
deduplicated archive with SHA-256
`b057954c12a969591a3d4451eb2a693b5f7ec6cbb611905a4c686733abedd422`.
#link("https://github.com/aburan28/crypto/actions/runs/39005111289")[GitHub Actions run 39005111289]
contains the two original provider artifacts.

== Frozen instance and replay

The exact curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis
$"GF"(2^83)$ with defining-polynomial integer
`9671406556917033397649557`, equation $y^2 + x y = x^3 + x^2 + 1$,
prime subgroup order `8569786107849059`, and cofactor `1128547018`.
The sixteen public points were sealed before this finalist panel;
`targets.jsonl` SHA-256 is
`407cd248797100214ff039193bb04bf89b5c2708a8cd339ca98f7aca98174062`.
The answer sidecar is opened only by the verifier.

The two provider ZIP SHA-256 values are
`50db2c8e06c884ad9f6a2731aa350dc8bde708d8ea0c52fddf740404874908a4`
and `737df73989d1a2e663ef3ea9253bb088b6aa92d2734070cae3a3d4b641b75c20`.
Local replay checked provider members, all factor-base orbits and subgroup
points, each four-point target relation and recovered scalar, exact source
and executed binary identity, event order, host telemetry and calibration.
The compressed archive stores all 3,802 logical file paths and uses four
distinct base-dump blobs for their 132 repetitions. A fresh reconstruction
passed every path hash in `raw_dedup_manifest.json`.

#pagebreak()

= Paired distribution and next experiment

#figure(image("finalist.svg", width: 100%), caption: [Sixteen independently verified public points per candidate. Thick colored bar: interquartile range; dark tick: median. The third panel divides each same-point rho online interval by the corresponding IC online interval.])

The candidate selector compared all 28 unordered pairs using exact two-sided
sign flips over $2^16$ assignments and a Bonferroni threshold of
$0.05 / 28 = 0.001785714$. Its smallest online $p$ is 0.045258. The
K1107/root11 versus K1107/root22 cold comparison has $p = 0.617950$.
The plotted data are pinned to the `validated_rows.json` and
`selection.json` SHA-256 values in `plot_data.json`; the editable graph
source is `build_visual.py`.

The next sealed ten-point experiment compares K1248/root11 and
K1280/root11 with K2543/root11, K3072/root11 and the predeclared
K1120/root11 fallback. It keeps same-point rho, two Tier 1 host strata,
alternating execution order and before/after rho calibration.
