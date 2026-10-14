# Root-finding decomposition engine in ecbench: results

**Decision: engineering.** The new oracles verify, replay identically, and
match the query stream of the existing table/subtract arms; they do not
change the asymptotic group-operation count of the decomposition, so
nothing here is an *advance*. Against Gröbner (`inherited-f4`) the root
finder's `S` is a complete count beside a lower bound and is not a
speedup. Against `mitm` / `subtract` the ops ratio is 1 (same stream);
wall time on this L0 macOS host is descriptive only and favours the root
finder. No security claim; no ECC2K-130 inference.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) (committed before any session ran) |
| engine ids | `pdp3-koblitz:m=2\|3,engine=rootfind`; `rootfind-subtract:m=2\|3` (see `docs/ecbench/README.md` §7) |
| specs | `docs/ecbench/specs/smoke-koblitz-rootfind.json`, `smoke-prime-rootfind.json`, `claim-koblitz-rootfind-n13.json`, `claim-koblitz-rootfind-n19.json`, `panel-prime-rootfind.json` |
| binary | `ecbench` release SHA-256 `1e972df4d8e769629ad4fd4a023214b05ee513a2cfea6f3ef51f2ca9e61fe531` |
| host | Apple M4 Pro, 14 logical CPUs, macOS 25.6.0 aarch64; env class `ECBENV2h59034d4cd61b`; `--cpus none --allow-busy` → every run **L0** |
| classification | **engineering** (same decomposition work, faster constant / cleaner charge); accounting vs Gröbner lower bounds |

## Complete sessions

| session | id | records | audit | replay cert (SHA-256 of audit JSON) |
|---|---|---:|---|---|
| [`sessions/smoke-koblitz`](sessions/smoke-koblitz) | `ECBS1hd86e5cfc2b77` | 42/42 verified | OK | `1cb597b686e68d783c8c370806686740ef316e862ebd191b6482470f2605c2b5` |
| [`sessions/smoke-prime`](sessions/smoke-prime) | `ECBS1hb8f2623c82e2` | 36/36 verified | OK | `bfb3f0207a2ae4f57ee8e176bd3cd05e44b8eb0d108340f70d4323e647188a1f` |
| [`sessions/claim-n13`](sessions/claim-n13) | `ECBS1hfd5bf1885227` | 224/224 verified | OK, 168/168 replays identical | `8e6a76596041d9a182e2584b9844ecf4846d827e32fbb1c9e3bd3d68e36a19ae` |
| [`sessions/panel-prime`](sessions/panel-prime) | `ECBS1h045961daac5c` | 768/768 verified | OK, 576/576 replays identical | `19d4aba584a99a48c46fd7b70262d31972a6d11770771826e2c5ac10a49d74ce` |
| [`sessions/claim-n19`](sessions/claim-n19) | `ECBS1h52e3269d5dab` | 224/224 verified | OK, 168/168 replays identical | `f3a3105be311cee69d7c2328a283a3c7632739136d832fd86d9de6e30998136d` |

Incomplete starts kept under [`abandoned/`](abandoned/) (SIGKILL / shell detach
at mid-session; status left as recorded).

## Boundary and unit

As in [`PROTOCOL.md`](PROTOCOL.md): floor `√(π/2A)` in `S`; reference
`rho.signed_frobenius_strong` (Koblitz) / `rho.negation` (prime); cold
`S` with every phase charged. Gröbner arms mark `lower bound` (word XORs
counted, not charged). Root-finding Koblitz arms are complete counts.

## Tables (mean `S`, L0)

Full `ecbench table` output: [`tables/`](tables/). Headline cells:

### Koblitz `E_0`, `n=13`, `d=5` (`claim-n13`, 8 public targets × 3 rounds)

| arm | oracle | mean S | S / strong-rho | lower bound |
|---|---|---:|---:|---|
| `rho-strong` | — | 14.118 | 1.000 | yes |
| `rho-frob` | — | 1.871 | 0.133 | yes |
| `ic-f4-m3` | `pdp3-koblitz:m=3,engine=inherited-f4,…` | 15.928 | 1.128 | yes |
| `ic-rootfind-m3` | `pdp3-koblitz:m=3,engine=rootfind` | 38.822 | 2.750 | **no** |
| `ic-f4-m2` | `pdp3-koblitz:m=2,engine=inherited-f4,…` | 16.642 | 1.179 | yes |
| `ic-rootfind-m2` | `pdp3-koblitz:m=2,engine=rootfind` | 28.383 | 2.010 | **no** |

`ic-rootfind-m3 / ic-f4-m3` ops = 2.44 [2.35, 2.53] is **diagnostic only**
(complete count over a lower bound). Descriptive solve-wall sums over the
24 measured runs: F4-m3 5.70 s vs rootfind-m3 0.13 s (~43×); F4-m2 0.45 s
vs rootfind-m2 0.12 s (~3.6×). Not admitted (L0).

Smoke Koblitz: `ic-rootfind-m2 / ic-mitm-m2` ops = **1.000** [1.000, 1.000]
(identical `S` and query stream).

### Koblitz `E_0`, `n=19`, `d=7` (`claim-n19`, 8 public targets × 3 rounds)

| arm | mean S | S / strong-rho | lower bound |
|---|---:|---:|---|
| `rho-strong` | 3.163 | 1.000 | yes |
| `ic-f4-m3` | 4.904 | 1.550 | yes |
| `ic-rootfind-m3` | 53.996 | 17.072 | yes |
| `ic-f4-m2` | 8.097 | 2.560 | yes |
| `ic-rootfind-m2` | 34.688 | 10.967 | yes |

Two sizes are not a scaling claim. Rootfind remains above rho and above the
Gröbner lower-bound `S` (complete charged table work vs uncharged solver XORs).

### Prime panel (18/20/22/24-bit `prime_search`, seeds of
`research/ecbench_all_candidates_20261003/spec-prime.json`)

| arm | 18-bit mean S (ref×) | 20-bit | 22-bit | 24-bit |
|---|---:|---:|---:|---:|
| `rho-neg` | 1.294 (1.00) | 1.344 | 1.480 | 1.080 |
| `ic-subtract` | 310 (239×) | 449 | 927 | 2775 |
| `ic-rootfind-m2` | 310 (239×) | ≈subtract | ≈subtract | ≈subtract |
| `ic-mitm` | 9.96 (7.7×) | 10.5 | 16.2 | 43.7 |
| `ic-rootfind-m3` | (see table) | | | |

`ic-rootfind-m2 / ic-subtract` ops = **1.000** [1.000, 1.000] over 96 pairs
(same candidates, same relations; GAE differs by one re-addition per
relation; `inversions` counted uncharged). Descriptive solve-wall:
subtract / rootfind-m2 ≈ 2.5× over the panel.

No IC arm is below rho on these curves.

## Hypotheses (from PROTOCOL)

- **H1:** every root-finding run verified; every audit replay `identical`. **Met.**
- **H2:** same stream as control oracles where the control completes. **Met**
  for `rootfind-m2` vs `mitm-m2` (Koblitz smoke) and vs `subtract` (prime panel).
- **H3:** root-finding solve wall below Gröbner at `m=3` on completed
  workloads. **Met descriptively** on claim-n13 (L0, not admitted).

## Inadmissible / caveats

- All wall figures are L0 on a shared macOS host (load often ≫ 1); do not
  cite them as admitted speedups.
- Two Koblitz sizes ≠ a scaling claim (AGENTS.md §5).
- Do not divide rootfind `S` by Gröbner `S`.
- Interrupted first attempts of claim-n13 are under `abandoned/`, not in the scoreboard.

## Reproduce

```bash
cargo build --release --bin ecbench
ECBENCH=target/release/ecbench bash research/ecbench_rootfind_engine_20261012/run.sh
```
