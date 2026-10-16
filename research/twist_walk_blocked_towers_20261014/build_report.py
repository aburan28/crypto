"""Lay out REPORT.typ from the recorded CSVs (layout only; computes nothing).

usage: python3 build_report.py catalog_towers.csv summary.csv class_ceiling.csv REPORT.typ
"""
import csv
import sys

catalog, summary, ceiling, out = sys.argv[1:5]
cat = list(csv.DictReader(open(catalog)))
summ = list(csv.DictReader(open(summary)))
null = [r for r in summ if r["kind"] == "null"]
cr = list(csv.DictReader(open(ceiling)))


def name(slug, n):
    return {"208": "c2pnb208w1", "368": "c2pnb368w1"}.get(n, f"{slug[:9]}…{slug[-8:]}")


def verdict(r):
    if r["nodes_d1"] == "0":
        return "no rational 3-isogeny"
    d = r["first_unblocked_depth"]
    if d == "none":
        return "no unblocked twist within 4 steps"
    return f"unblocked at depth {d}"


rows = []
for r in cat:
    rows.append(
        f"[{name(r['sample'], r['N'])}], [{r['N']}], [{r['l']}], [{r['n']}], [{r['start_m']}], "
        f"[{r['nodes_d1']}, {r['nodes_d2']}, {r['nodes_d3']}, {r['nodes_d4']}], [{verdict(r)}], "
        f"[{r['min_m_unblocked_le4'] or '—'}]"
    )
cat_table = ",\n  ".join(rows)

nrows = []
for r in null:
    nb = r["has_neighbour"]
    esc = r["escape_le4"]
    gen = r["escape_with_generic_m"]
    pct = lambda s: f"{100 * int(s.split('/')[0]) / int(s.split('/')[1]):.1f}%"  # noqa: E731
    nrows.append(
        f"[{r['N']}], [{r['l']}], [{r['n']}], [{r['sample']}], [{pct(nb)} ({float(r['has_neighbour_ci_lo']):.3f}–{float(r['has_neighbour_ci_hi']):.3f})], "
        f"[{pct(esc)} ({float(r['escape_le4_ci_lo']):.4f}–{float(r['escape_le4_ci_hi']):.4f})], [{gen}]"
    )
null_table = ",\n  ".join(nrows)

mixed_l1 = [r for r in cr if r["l"] == "1"]
mixed_l1_total = sum(int(r["classes_mixed"]) for r in mixed_l1)
lo_frac = min(int(r["curves_blocked_in_mixed_classes"]) / int(r["curves_blocked"]) for r in cr if r["l"] != "1")
hi_frac = max(int(r["curves_blocked_in_mixed_classes"]) / int(r["curves_blocked"]) for r in cr if r["l"] != "1")
nl1 = len(mixed_l1)

typ = f'''#set page(paper: "a4", margin: (x: 1.8cm, y: 1.5cm))
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set block(above: 1em, below: 0.5em)

= Leaving a blocked twist tower: 3-isogeny walks and the class ceiling

_2026-10-14. Frozen protocol `PROTOCOL.md` (sha256 `{open("manifest.json").read().split('"protocol_sha256": "')[1][:16]}…`). Native instrument `examples/twist_walk_blocked_towers.rs`. Both positive controls passed before any data. Scope: 5 catalog curves, 13 towers, depth 4, degree 3 only._

== Catalog result (exact computation)

Twist of each of the 13 blocked towers; depth counts are 3-isogeny steps, Frobenius free. `m` is the MMT magic number of the start; `m'` is the least magic number among unblocked nodes within depth 4.

#table(
  columns: (auto, auto, auto, auto, auto, auto, 1fr, auto),
  align: (left, right, right, right, right, left, left, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*Curve*], [*N*], [*l*], [*n*], [*m*], [*nodes d1,d2,d3,d4*], [*verdict*], [*m'*],
  {cat_table},
)

Of 13 towers, 3 are unblocked after one step (`N=24, l=2`; `c2pnb368w1, l=2` and `l=46`). 7 have no rational 3-isogeny at all (three curves, including `c2pnb208w1`). 3 have neighbours and stay blocked: the two `l = 1` towers and `c2pnb368w1, l = 23`.

== Random blocked starts, for context (not tested against the catalog)

2000 random starts per configuration, same field, `Tr_(K\\/k)(b) = 0` forced; 95% Clopper–Pearson intervals. The last column counts escapes whose least `m'` equals `n` (the generic maximum).

#table(
  columns: (auto, auto, auto, auto, 1fr, 1fr, auto),
  align: (right, right, right, right, left, left, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*N*], [*l*], [*n*], [*starts*], [*has a 3-isogeny*], [*escaped within 4 steps*], [*escapes with m' = n*],
  {null_table},
)

#figure(image("escape_by_config.svg", width: 88%), caption: [Escape rate of random blocked starts within 4 steps (blue, 95% interval) and the catalog towers (orange). Source: `summary.csv`.])

== Class ceiling (exact computation, finite)

The twist-side isogeny class of `b` is `{{b' : t(b') = t(b)}}`. For every field `N = 8…22` and every `l | N` with `n = N/l` even we counted the trace classes holding a blocked curve and those that also hold an unblocked one (`class_ceiling.csv`). For `l = 1` no class is mixed in any of the {nl1} field sizes ({mixed_l1_total} mixed classes in total), so no isogeny of any degree leaves an `l = 1` blocked tower in those fields. For `l >= 2`, {100 * lo_frac:.0f}% to {100 * hi_frac:.1f}% of blocked curves lie in a mixed class, so an escape exists inside the class.

#figure(image("class_ceiling.svg", width: 88%), caption: [Share of blocked curves whose isogeny class also contains an unblocked curve, by `l`. Source: `class_ceiling.csv`.])

== Claims and scope

- *Exact computation:* the 13 catalog outcomes above; the class counts for `N = 8…22`.
- *Observed on a sample:* 0 of 2000 random starts escape an `l = 1` tower in each of four fields (upper limit 0.0018), as the class counts require; 36–40% escape at `l = 2`, 65% at `N = 18, l = 3`, and 74–75% at `l >= 13`.
- *Not tested:* that the catalog differs from random starts. Five independent graphs cannot support that test at the corrected level (`PROTOCOL.md`). One catalog escape is far from generic: `c2pnb368w1` at `l = 2` reaches `m' = 8` where random escapes all have `m' = n = 184`; the other two escapes have `m' = n`. That contrast is descriptive.
- *Conjecture, not proved here:* `Tr_(K\\/F_2)(b)` is constant on each trace class for all even `N`; the exhaustive counts support it for `N <= 22`.
- *Not shown:* that any twist is attackable. An unblocked twist tower removes the Lemma 6 obstruction only; the twist's group order and cofactor, the cover's genus and the descent cost are separate. Degrees other than 3 are not walked, and depth is capped at 4.
'''
open(out, "w").write(typ)
print("wrote", out)
