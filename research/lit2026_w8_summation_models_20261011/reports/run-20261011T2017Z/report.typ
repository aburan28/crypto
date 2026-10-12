#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#set par(justify: false)

= Run report: lit2026 W8 diagnostics, run-20261011T2017Z

*Started* 2026-10-11 20:17 UTC. *Host* shared Apple M4 Pro, macOS (Darwin 25.6.0), load average 72–360
during the parts (`uptime` per round in every JSON). *Isolation* L0. *Harness*
`examples/lit2026_w8_summation_models.rs` on branch `lit2026/w8-summation-models`; library byte-identical
to program baseline `4b6e09b78` (`git diff 4b6e09b78 HEAD -- src Cargo.toml Cargo.lock` empty).
This report describes what ran; interpretation is in the study README.

== Parts, commands, binaries

Every part: `--rounds 5 --targets 8`, arms interleaved with the order rotating per round, seed `0x5EED`.
Binaries are successive builds of the same example; SHA-256 of each and what changed:

#table(columns: 4, inset: 4pt,
  [*file*], [*binary (sha256 prefix)*], [*env*], [*status*],
  [`cert.json`], [`420148c1` (re-issued by `da887b81`, byte-identical result)], [—], [complete],
  [`lucas.json`], [`420148c1`], [default], [complete (n = 9, 12, 15)],
  [`lucas-l34.json`], [`b42e8494`], [default], [complete (n = 9, 12); word ops equal `lucas.json`],
  [`wu-first-binary-failed.json`], [`420148c1`], [default], [failed: binary Edwards fit used 3 points that repeat ±P (bug, fixed)],
  [`wu.json`], [`da887b81`], [default], [complete],
  [`eg-aborted-first-binary.stderr`], [`420148c1`], [default], [stopped (rerun with ladder-only env)],
  [`eg.json.partial`], [`5216edd7`], [`KIC_F4_MAX_DEGREE_ONLY=0`], [9 of 10 cells; K1/F2^17 m=3 stopped after 2 h (SIGTERM)],
  [`frob-aborted-first-binary.stderr`], [`420148c1`], [default], [stopped: aux arm at 28 vars on the top-degree-only path, > 17 CPU-min in cell 1],
  [`frob-aborted-w8bin3.*`], [`5216edd7`], [`=0`], [stopped to add the engine flag; 3 cells recorded],
  [`frob.json`], [`99a4a31b`], [`=0`], [complete, MatrixF4],
  [`frob-inherited.json`], [`99a4a31b`], [`=0`], [complete, `--engine inherited`],
)

Builds: `420148c1` first; `da887b81` adds the binary-Edwards fit against `X = x + sqrt(b)/x` (Wu part only);
`5216edd7` adds per-cell checkpoints, the env record and the reduced `frob` cell list; `99a4a31b` adds
`--engine` (default unchanged) and `resolution_degrees_built`; `b42e8494` adds `--arms`, `--no-degrees`,
root digests and `compare`. No build changed the computation of a part it did not target.
Driver logs: `run.log`, `run_eg3.log`, `run_wu2.log`, `run_frob4.log`.

== Workloads

Koblitz curves `y^2 + xy = x^3 + a x^2 + 1` from `KoblitzCurve::new(a, n)`: K0/F2^9 `icv1-f2m9-t5-81e744be`
(r = 127), K1/F2^15 `icv1-f2m15-t275-b7f03703` (r = 211), K0/F2^15 `icv1-f2m15-tm275-2d22ff5d` (r = 751),
K1/F2^17 `icv1-f2m17-tm101-00378d4e` (r = 65 587), K0/F2^23 `icv1-f2m23-t5197-69e76b73` (r = 2 095 853),
K1/F2^23 `icv1-f2m23-tm5197-1f85e9e1` (r = 4 196 903); targets uniform in the prime-order subgroup.
Wu part also `icv1-f2m15-t5-d9fb183e` and `icv1-f2m17-t209-0eba1a1e` (a = 0, b != 1), targets uniform affine
points. Factor bases: `eg`/`frob` use the divisor-of-`x^n - 1` subspace of `koblitz_symmetrised`
(dimension ell in the JSON); `wu` uses `<1, z, ..., z^(ell-1)>`; `lucas` uses `<1, ..., z^(l-1)>` with
`n = 3l` and random target abscissae.

== Recorded values

#figure(image("eg_arms.svg", width: 100%), caption: [E/G check, `eg.json.partial`: median word XORs per
solve and Boolean unknowns per arm and cell (8 targets x 5 rounds; MatrixF4 cap 4; `KIC_F4_MAX_DEGREE_ONLY=0`).])

#figure(image("wu_models.svg", width: 100%), caption: [Wu part, `wu.json`: S3 terms and total degree per
coordinate, and median word XORs per curve (8 targets x 5 rounds; MatrixF4 cap 4).])

#figure(image("frob_lucas.svg", width: 100%), caption: [`frob.json`, `frob-inherited.json`, `lucas.json`:
median word XORs per arm; root sets identical between arms in every cell.])

== Gates recorded

- `eg`: `w`-root sets of `sym` (parity dropped), `eg` and `u-plain` (image under AS) identical on all 72
  targets; `sym` verdicts equal exhaustive enumeration on all 72; root sets stable across rounds;
  `b''` of `E/<T>` read and checked on 32 relations = 1 on every cell.
- `wu`: every model polynomial vanished on 32/32 random relations; every descended root set equal to
  brute force over `V x V` on 32/32 targets per model; Wu equation on 64/64 sampled points per curve;
  binary Edwards curve equation on 64/64 points and `w` invariant under `T` on 16/16.
- `frob`: projected root sets identical on 32/32 targets, both engines.
- `lucas`: projected root sets equal brute force over `V^3` on 24/24 targets, both splits.
- `cert`: certificates IDC1hc7b594938b50f160 (S3) and IDC1hcf83176ae9a93715 (S4) issue and check (32 points each).

Missing data: K1/F2^17 m=3 in the `eg` part (stopped). No other requested cell is missing.
