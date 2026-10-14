#set document(title: "MMT Lemma 6 in the Jacobian-transfer screen")
#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(size: 10pt)
#set par(justify: true)
#show link: set text(fill: rgb("#1c5cab"))
#show table: set par(justify: false)
#show table: set text(size: 9pt)

#text(size: 15pt, weight: "bold")[MMT Lemma 6 in the Jacobian-transfer screen]
#linebreak()
#text(fill: rgb("#6b6a63"))[2026-10-13 · `research/jacobian_transfer_lemma6_20261013/` · PR \#1772 (aburan28/crypto)]

== Question

Take a binary curve $E: y^2 + x y = x^3 + a x^2 + b$ over $K = bb(F)_(2^N)$ with $"Tr"(a) = 1$ and an even relative degree $n = N slash l$. What are the GHS magic number and the descent route over $k = bb(F)_(2^l)$? Before this change the screen (`src/cryptanalysis/jacobian_transfer/`) did not answer that question. It priced two heuristic rows, $m$ and $m+1$, and marked the curve `unresolved` when they disagreed.

== Answer (derived; exact on each curve)

- $a$ never changes the magic number $m = dim "Span"{(1, sigma^i (sqrt(b)))}$ (Maurer–Menezes–Teske, LMS J. Comput. Math. 5 (2002), eq. (5)).
- $a$ decides whether the descent exists. MMT Lemma 6, condition (8): $n$ odd, or $"Tr"(a) = 0$, or $"Tr"_(K slash k)(b) != 0$.
- Why: when $n$ is even and $"Tr"_(K slash k)(b) = 0$, the sum of the $n$ Artin–Schreier conjugates of $z^2 + z = x + a + b slash x^2$ reduces to the constant $"Tr"_(K slash k)(a)$. Its $℘$-root lies in $bb(F)_(q^2)$, and lies in $k$ only when $"Tr"(a) = 0$. Otherwise every lift of the $q$-Frobenius has order $2n$, and there is no descent.
- The screen now gives one exact row per tower; a row where (8) fails is `descent_obstructed`.

The full derivation and the list of changed code are in `docs/audits/GHS_MAGIC_NUMBER_MMT_CORRECTION_20261013.md`.

== What changed

Sources:
- *Before:* `docs/curves/jacobian-transfers.json` at commit `b8ce3f697`.
- *After:* the same file at this report's commit.
- Registry SHA-256 `9e420a7f…1219e872`, `--budget 200000`.

The values are copied by `extract_rows.py`. Magic numbers, genera and condition (8) are exact invariants. The ratio $"log"_2("rho") - "log"_2("route")$ is a formula under the frozen cost model, not a measurement.

#table(
  columns: (auto, auto, 1.6fr, 1.4fr, 1.3fr),
  stroke: 0.4pt + rgb("#c3c2b7"),
  inset: 4pt,
  align: left,
  table.header([*curve*], [*$r$, bits*], [*GF(4) row before*], [*GF(4) row after*], [*verdict before → after*]),
  [c2pnb176w1], [160.0], [$m=8$ lower bound, heuristic, +23.64; $m=9$ upper bound, −5.76], [$m=8$ exact, genus 128, (8) holds, +23.64], [unresolved → formula_advantage],
  [c2pnb272w1], [256.0], [$m=8$ lower bound, heuristic, +72.20; $m=9$ upper bound, +41.32], [$m=8$ exact, genus 128, (8) holds, +72.20], [formula_advantage → formula_advantage],
  [c2pnb304w1], [288.0], [$m=8$ lower bound, not carried; $m=9$ upper bound, +57.56], [$m=8$ exact, genus 128, not carried (256 < 288 bits)], [unresolved → formula_marginal (GF($2^76$), $m=4$, genus 8, +2.00)],
)

#block[#set par(justify: false)
#text(size: 8.5pt)[ICV1 slugs:
- c2pnb176w1: `icv1-f2m176-t38161931374797611715595947-5639a6d4`
- c2pnb272w1: `icv1-f2m272-t114182707820475866940422329614334265958971-985955bd`
- c2pnb304w1: `icv1-f2m304-tm1949050066876782037155448244017940336576672053-02fbe089`]]

Across all 2091 registry curves, these three are the only verdict or GF(4) row changes.

The isogeny-extension rows of these three curves are gone. For $"Tr"(a) = 1$ and even $n$ they now count only targets that satisfy (8), and none remain below the direct $m$. Every removed row was `not_carried` or `boundary`.

Condition (8) fails on one registry curve, `icv1-f2m18-t999-40751283` ($b in bb(F)_2$, $"Tr"(a) = 1$), at $l = 9, 3, 1$. Its rows were already not carried, and its verdict stays `boundary`.

#figure(
  image("gf4_rows.svg", width: 100%),
  caption: [GHS rows over GF($2^2$) before (`b8ce3f697`) and after, plus c2pnb304w1's new best route (`gf4_rows.svg` from `gf4_rows.csv`). Hollow markers are not carried; `*` marks each curve's best route. The dashed line is the 8-bit advantage margin.],
)

#figure(
  image("binary_curves.svg", width: 100%),
  caption: [Best carried route of every binary-field curve in the file (160 curves): 155 `boundary`, 3 `formula_marginal`, 2 `formula_advantage` (`binary_curves.svg` from `binary_curves.csv`). The two curves whose verdict this change moved are labeled.],
)

== Checks run

- `cargo test --release --lib`: 4065 passed, 0 failed, 114 ignored.
- `mmt_condition_8_matches_the_relative_trace_of_b`: exhaustive over $bb(F)_(2^12)$ on all six towers. It compares condition (8) with $sum b^(q^i)$ computed directly and checks the kernel size $q^(n-1) - 1$. It also checks that $m$ and the type do not depend on $a$, and that `ec_trapdoor::mmt_condition_8` agrees.
- `ec_trapdoor` tests: `mmt_condition_8_counts_n12`, \ `audit_curve_puts_descending_towers_first`, \ `trapdoor_curves_satisfy_condition_8`.
- `jacobian_transfers check` and the `JACOBIAN_TRANSFERS.md` table: current.

== Scope and limits

- These are formula-level prices on exact invariants. No cover is constructed, no map's action on the subgroup is verified, no relation is collected, and nothing is timed.
- No speedup, key-recovery or deployed-curve security claim is made. The `formula_advantage` rows are leads with every obligation in `docs/endomorphism-rules.md` open.
- Condition (8) is used as MMT state it. The exhaustive test checks this implementation against the definition; it does not re-prove the lemma.

== Graphs checked

- The 2026-10-10 scan's `chart.svg` and `report.pdf` are a dated 321-curve snapshot. They are left unchanged, and its README has a dated correction note pointing here.
- The registry-wide `jacobian_transfers chart` labels every advantage and marginal curve (846 + 47 of them toy extension-field curves), and those labels overlap at 2091 curves. Figure 2 is drawn from the same file instead.
