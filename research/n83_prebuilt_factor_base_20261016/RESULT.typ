#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)

= N83 prebuilt factor-base result

*Protocol:* `PROTOCOL.md`, frozen at commit
`2d3a7f3723925311e74a29ff1cdc7c74a0a97f50`  \
*Curve:* `icv1-f2m83-tm6151469093347-debefd74` primary;
the degree-83 `a=1` curve is a separate diagnostic.  \
*Host:* shared macOS arm64, L0, 2026-10-16.

#image("reports/run-20261016T2352Z/panel-times.svg", width: 100%)

Three fresh native builds each completed all 54 declared objects. Each
panel contains 83,328,324 compressed bytes. The canonical listing of point
sets, plain bytes and compressed bytes is identical across the three builds,
one later reconstruction and the historical content-addressed S3 panel.
Its SHA-256 is
`e0079cf8ad843084243a1319de1aefa7ccf65ed1f71a92aaddd588761e514e3d`.

#table(
  columns: (auto, auto, auto, auto),
  inset: 5pt,
  align: (right, right, right, right),
  [*Panel*], [*Build, internal s*], [*Outer real s*], [*Generic replay s*],
  [1], [59.518], [60.07], [46.66],
  [2], [74.022], [74.73], [52.15],
  [3], [169.555], [172.78], [31.43],
)

The three generic multi-limb replays passed 54 checks each, using the
original manifests and byte-identical reconstructed objects. The first
three replay attempts failed because temporary object paths had been
cleared; those failure receipts remain in `runs/`. The timing wrapper
exited 1 after each completed build because sandboxed `sysctl kern.clockrate`
was unavailable; peak RSS is missing. The 2.85-fold build interval spread
prevents an isolated wall-performance inference. Build and replay perform
different work and their times are not a speed ratio.

For the primary K=600 orbit base, the exact multiset counting ceilings
for a uniform nonidentity subgroup target are 2.051462e-15, 6.810989e-11
and 1.695987e-6 at two, three and four summands. They bound one-query
full-smooth probability under that model; actual relation yield and useful
rank remain unmeasured. The stored dimension-12 structured base has 2,001
signed-Frobenius columns and an exact four-summand count/order ratio of
0.000209791390, also a ceiling rather than measured yield.

*Decision:* Ahead-of-time construction and checked reuse are demonstrated
at the factor-base stage. Selecting the best symmetric base at two to four
summands requires a matched natural-query oracle and rank comparison of
the near-size-matched K=2,048 explicit-orbit and 2,001-column structured
bases, followed by complete cold accounting and controls. See `RESULT.md`,
`REVIEW.md`, `VALIDATION.md`, `runs/` and
`problems/n83-prebuilt-base-quality-m2-m4.md`. No canonical performance
graph changes: there is no admitted end-to-end improvement measurement.
