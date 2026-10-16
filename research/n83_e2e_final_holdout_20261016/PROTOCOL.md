# N83 finalist holdout: sixteen frozen public points

This fixture freezes sixteen new N83 subgroup points before the adjacent-K
and precision results are combined with the completed cross-layout screen
into a finalist choice. The curve ID is `EC1N83Ckb1a1haacb69ee9f37`:
the field is GF(2^83) in a polynomial basis with defining-polynomial integer
`9671406556917033397649557`, the curve is `y^2+x*y=x^3+x^2+1`, and the
prime subgroup order is `8569786107849059`. The scalar rule is
`1 + little_endian(SHA256(domain || u32_le(index))) mod (r - 1)` for indices
0 through 15. The public points are the measured inputs; the answer sidecar
is used only for independent scalar verification.

`generate.py` pins the precision fixture generator and both precision input
file hashes. It replays that generator's earlier target exclusions, checks
the sixteen precision points, and rejects any repeated, off-curve or
out-of-subgroup point before sealing this new target file. The committed
`fixture_manifest.json` records the exact public-point and sidecar hashes.
The public-point SHA-256 is
`407cd248797100214ff039193bb04bf89b5c2708a8cd339ca98f7aca98174062`;
the answer-sidecar SHA-256 is
`fe771284a3b0e5c0301309f678751f0630f0fa23844dd0ce1dbc9fa47f243076`.
Run `python3 research/n83_e2e_final_holdout_20261016/generate.py
--verify-sealed` before using the fixture in a measured job.

The next frozen protocol will name the exact finalist candidate manifests
after the pending screens complete. It must carry two separate comparisons:
the primary single-target online interval after reusable preparation, paired
with same-point rho, and the supplementary complete one-target cold process
interval. K3072/root11 is the current independently replayed Tier 1 online
reference; the K1104/K1120 root11 and root22 families are the current cold
finalists. The exact executed-source candidate IDs will be frozen before
running this fixture. Each comparison will use fresh processes on these same
public points, symmetric A/B controls, separately timed phases, and
alternating same-point rho solves within a Tier 1 resource envelope. It will
fix paired advancement rules for both metrics before any timings on these
points, preserving inconclusive selections as null.
