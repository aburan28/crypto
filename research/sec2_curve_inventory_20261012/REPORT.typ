#set page(paper: "a4", margin: 17mm)
#set text(font: "Arial", size: 9pt)
#set par(justify: true)
= SEC 2 curves: identity and structural traits
2026-10-12 · primary-standard crosscheck · 33 named curves

The 2000 edition lists 33 curves. The 2010 edition retains 20 of them; 13 are historical entries. Every field, coefficient, uncompressed base point, subgroup order, and cofactor in both editions matched the exact crypto registry record.

#figure(image("sizes.svg", width: 100%), caption: [Field size and edition membership. Circle: prime field; square: binary field. Blue: listed in 2010. Amber: 2000 only.])

#pagebreak()
== Curve index
Each published base point has order n, and the full curve has order h × n. Full exact identities, generators, n values, seeds, small-prime splitting, conductors, embedding and twist records, cover and transfer screening, and all status fields are in `inventory.json`. F in the last column means the SEC 1 advisory bound fails; ? means the bounded n±1 check was inconclusive.

#set text(size: 7.3pt)
#table(columns: (1.3fr, 0.8fr, 0.8fr, 0.5fr, 0.9fr, 0.8fr, 0.9fr), inset: 3pt, stroke: 0.35pt + rgb("dddddd"),
[Curve], [Edition], [Field], [h], [CM status], [Twist], [n±1],
  [secp112r1], [2000 only], [Fp/112], [1], [probable], [proved], [F/?],
  [secp112r2], [2000 only], [Fp/112], [4], [probable], [proved], [F/?],
  [secp128r1], [2000 only], [Fp/128], [1], [bounded], [proved], [F/F],
  [secp128r2], [2000 only], [Fp/128], [4], [probable], [probable], [F/F],
  [secp160k1], [2000 only], [Fp/160], [1], [-3 (proved)], [probable], [F/F],
  [secp160r1], [2000 only], [Fp/160], [1], [probable], [probable], [?/F],
  [secp160r2], [2000 only], [Fp/160], [1], [probable], [probable], [F/F],
  [secp192k1], [2010], [Fp/192], [1], [-3 (proved)], [probable], [F/F],
  [secp192r1], [2010], [Fp/192], [1], [probable], [unknown], [F/?],
  [secp224k1], [2010], [Fp/224], [1], [-3 (proved)], [unknown], [?/F],
  [secp224r1], [2010], [Fp/224], [1], [probable], [unknown], [F/?],
  [secp256k1], [2010], [Fp/256], [1], [-3 (probable)], [probable], [F/F],
  [secp256r1], [2010], [Fp/256], [1], [bounded], [probable], [F/F],
  [secp384r1], [2010], [Fp/384], [1], [bounded], [probable], [?/F],
  [secp521r1], [2010], [Fp/521], [1], [bounded], [probable], [?/?],
  [sect113r1], [2000 only], [F2m/113], [2], [proved], [proved], [F/F],
  [sect113r2], [2000 only], [F2m/113], [2], [probable], [proved], [F/F],
  [sect131r1], [2000 only], [F2m/131], [2], [probable], [probable], [F/?],
  [sect131r2], [2000 only], [F2m/131], [2], [probable], [probable], [?/F],
  [sect163k1], [2010], [F2m/163], [2], [-7 (proved)], [unknown], [F/F],
  [sect163r1], [2010], [F2m/163], [2], [bounded], [unknown], [F/F],
  [sect163r2], [2010], [F2m/163], [2], [bounded], [probable], [F/?],
  [sect193r1], [2000 only], [F2m/193], [2], [probable], [unknown], [F/F],
  [sect193r2], [2000 only], [F2m/193], [2], [probable], [unknown], [?/?],
  [sect233k1], [2010], [F2m/233], [4], [-7 (proved)], [unknown], [F/F],
  [sect233r1], [2010], [F2m/233], [2], [bounded], [unknown], [F/?],
  [sect239k1], [2010], [F2m/239], [4], [-7 (proved)], [unknown], [F/?],
  [sect283k1], [2010], [F2m/283], [4], [-7 (probable)], [probable], [F/?],
  [sect283r1], [2010], [F2m/283], [2], [bounded], [unknown], [F/?],
  [sect409k1], [2010], [F2m/409], [4], [-7 (proved)], [unknown], [F/F],
  [sect409r1], [2010], [F2m/409], [2], [bounded], [unknown], [?/?],
  [sect571k1], [2010], [F2m/571], [4], [-7 (probable)], [unknown], [F/F],
  [sect571r1], [2010], [F2m/571], [2], [bounded], [unknown], [?/?],
)

== Sources, method, and status
SEC 2 v1.0 (2000): `https://www.secg.org/SEC2-Ver-1.0.pdf`
SHA-256 `d1b16728ad83888fd656d16b99dc71bcd5541d42d848ffd0de7c62c19010d8c3`.
SEC 2 v2.0 (2010): `https://www.secg.org/sec2-v2.pdf`
SHA-256 `87b8f3703364ed5b21ba8582e411cc0cbf477bcaa3f4f45e0d6580d1c00d9952`.

The two PDFs supply 33 and 20 exact named parameter blocks. `build.py` compares field prime or reduction polynomial, coefficients, uncompressed base point, prime subgroup order, and cofactor with the registry. All 53 blocks match. Seeds are transcribed but the ANSI seed-to-curve derivation is not independently replayed.

Existing `curve_traits` records give subgroup, CM/conductor, endomorphism bounds, small-prime splitting, subfield, embedding, and twist data. Cover and Jacobian-transfer screening rows match all 33 exact identities. Every cover has subgroup transport `not_tested`; the transfer formulas are screening models.

The bounded n±1 scan certifies failure of the later SEC 1 19/20-factor advisory for 25 n−1 values and 19 n+1 values. Among curves retained in 2010, 15 n−1 and 10 n+1 values also fail, so failure of this advisory condition alone cannot account for the 193-bit omission. The other cases remain unresolved. The trait and transfer files have older whole-registry pins, but all 33 selected rows match current exact identities and order/trace.

== Evidence boundaries
This is a primary-parameter check and a structured view of existing catalog traits. A `probable` or `bounded` trait is not promoted to proved. The n±1 scan uses trial division through 100,000 and proves failure of the SEC 1 advisory threshold only when the residual is already below it. Covers and transfer screening do not verify subgroup transport or measure attack cost. Full status-bearing records are in `inventory.json`; source PDF copies are in `sources/`.
