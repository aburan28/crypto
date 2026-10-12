#set page(paper: "a4", margin: 1.8cm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#text(13pt, weight: "bold")[Run report: run-20261012T0002Z-binary-validate]

*What ran.* `./target/release/examples/binary_gls_validate --m 31 --count 3 --hmax 320 --lmax 61 --points 24 --seed 20261012 --json research/isogeny_loop_glv_20261012/binary_validate.json`. Final source (clippy-clean). rustc 1.97.0, release, cloud container shared with another agent's builds, not isolated. Start: `../../binary_validate_start_utc.txt`. Wall time: `../../binary_validate.log`.

*Method.* Curves over F#sub[2#super[31]] are built from H#sub[D] mod 2, with H#sub[D] computed natively; the field is GF(2#super[62]). The walk takes an odd ℓ-isogeny from Kohel kernel polynomials, keeps it only if its codomain is on the crater, and closes it with F#super[kk]. The GLS twist is a#sub[2] + δ with Tr(δ) = 1. Checks run x-only on 24 points of prime order r on E′:
- x(ψP) = x#super[q] against x([μ]P);
- (R(x))#super[2#super[kk]] against x([λ(±s)]P).

*Runs requested / completed.* Three attempts:
+ hmax 160, ℓ ≤ 41: found 1 discriminant before |D| reached 10#super[7]. PASS. Kept as `../../binary_validate_run1_hmax160.*`.
+ hmax 320, ℓ ≤ 61: found 2, both PASS, 11 min 10 s. Its output stops after the two PASS rows: no verdict line and no JSON. The cause was not determined; the run was under `timeout 900` inside a shell that the tool moved to the background. Its printed rows are quoted in `../../binary_validate_cmd.txt`.
+ The re-run recorded here.

*Recorded values* (`../../binary_validate.md`):
#table(columns: 6,
  [*D*], [*h*], [*loop*], [*H#sub[D] bits / worst err*], [*ψ*], [*loop acts as*],
  [−39583], [62], [ℓ = 31, 11 F], [4162 / 2#super[−200]], [24/24], [±λ(s) 24/24],
  [−482623], [186], [ℓ = 31, 12 F], [14081 / 2#super[−473]], [24/24], [±λ(s) 24/24],
)

#figure(image("../../figures/status.svg", width: 92%), caption: [Where this check sits (amber, bottom middle). Diagram, not data: two pass/fail rows give no comparable quantities to plot.])

*Observations only.* Interpretation: `../../README.md` §4.
