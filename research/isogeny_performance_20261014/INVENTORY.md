# Implemented isogeny algorithms and measurement boundaries

The repository has considerably more than six methods. This inventory names
families without claiming they solve the same problem. Module paths below are
relative to `isogeny_algos/src/`; the broader `src/isogeny/` and cryptanalysis
walkers are consumers or older, separately bounded implementations.

| Input -> output | Implementations | Existing native `bench` groups | New exact-output series |
|---|---|---|---|
| Kernel points -> normalized map/codomain | Vélu cyclic/general, x-only Vélu, Montgomery affine/projective Vélu, Weierstrass and Montgomery square-root Vélu (`kernel/`) | `kernel`, `kernel2`, `p1kernel`, `big` | Next constructor series; preserve point-generation and map-evaluation boundaries |
| Kernel polynomial -> normalized map | Kohel odd/even (`kernel/kohel.rs`), binary Kohel | `kernel`, `kernel2`, `char2`, `big` | Independent reference in first series; never compare its supplied-kernel input directly to a finder |
| Isogenous curve pair -> normalized kernel/map | Linear algebra, Stark 1972, Atkin 1992, Atkin modular composition, Elkies 1992, Elkies 1998, fastElkies, fastElkies prime (`find/bmss.rs`); older Padé implementation (`find/elkies.rs`) | `bmss`, `big` | All eight BMSS methods at 192/224/256/511 bits; exact kernel, numerator, denominator and codomain fingerprints |
| Curve and prime degree -> rational isogenies | Division-polynomial factoring (`find/divpoly.rs`); modular-polynomial roots plus Elkies/BMSS (`find/elkies.rs`) | `find`, `char2`, `bmss` | Next finder series; charge modular-polynomial creation separately and include it in cold totals |
| Prime degree -> modular polynomial | Dense linear algebra, Hecke/Newton q-expansion and integer CRT (`find/modpoly.rs`); volcano/CRT Sutherland (`find/sutherland.rs`) | `phi`, `nc_suth` | Next polynomial series; identical coefficient fingerprints and setup work |
| Curve -> point count | Schoof–Elkies–Atkin, isogeny cycles, Atkin extension-field eigenvalues (`find/sea.rs`) | `sea`, `nc_atkin` | Auxiliary finding cost; keep point count distinct from constructed maps |
| Curve pair -> prime-degree edge path | Galbraith BFS, GHS, Kohel volcano, Couveignes MITM, Delfs–Galbraith (`path/`) | `path` | Next path series; exact endpoint and edge receipts, failed searches retained |
| Class-group action/path | Galbraith–Stolbunov weighted walk, CSIDH stepwise/batched/tree, ideal-order/relation search (`path/csidh.rs`, `relation.rs`) | `csidh`, `v2path`, `relation` | Next action series; complete actions and setup counted, planted endpoints disclosed |
| Special models / characteristic | Edwards/Huff, Hessian, radical 3/5/7, binary fields, characteristic-3 general Weierstrass | `models`, `radical`, `char2` and dedicated tests | Separate model-specific series; no unsupported adapter implied |
| Composite-degree chain | Naive/balanced/cost-model odd chains; Montgomery 2/4-isogeny power-of-two chains | `chain`, `twopow` | Next chain series; same points pushed and strategy-generation boundary |
| Quaternion ideals / supersingular paths | KLPT, Brandt/class sets, Deuring (`quat/`) | `quat` | Separate ideal and map boundaries; small-characteristic Deuring limits retained |
| Abelian varieties | Richelot, gluing/splitting, dimension-2 and general theta/Kani chains | `genus2`, `theta`, `theta4` | Separate genus/dimension series; not elliptic prime-degree constructors |
| Named NIST curve and large prime -> one descended eigenline map | Extension-field torsion, Frobenius eigenspace and descended kernel in `tools/isogeny-cli/src/constructor.rs`; exact independent replay in `replay.rs` and `verification_*` | Dedicated `isogeny search` / `walk` | Required whole-search series on published P192/P224 at 1117/2341 and larger degrees; construction and replay charged separately |
| CM curve / ideal loop | Hilbert class polynomial and curve construction (`cm.rs`); closing CM isogeny loops (`gls_loop.rs`) | Dedicated module tests / examples | Required construction and loop series; no kernel-only or endomorphism-existence timing substituted |

The root crate also contains earlier implementations in `src/isogeny/`: explicit
kernel and odd-degree/degree-2 Vélu, naive dual/composition, CM trace/discriminant,
binary quadratic forms/class groups, graph search and volcano traversal. Root
cryptanalysis consumers include binary isogenies and the P256 walkers/covers.
They need separately pinned native entrypoints and fixtures rather than being
silently equated with the standalone crate's methods. The root test suite
checks their correctness but is not their performance baseline.

The broad existing harness exposes 24 named groups (plus aggregates). Its
historic results are not new CPU baselines. Every subsequent series will use
fixed fixtures and hashes with its own correctness boundary; unsupported or
resource-censored rows remain explicit. The first series implements BMSS
counting and output identity. Other rows are inventoried, not measured by
that series, and remain in the campaign backlog. The Linux matrix also invokes
all 24 explicit legacy benchmark groups under the native isolation wrapper.
Those exploratory rows retain failed verification flags and lack full output
digests; they cannot establish an optimization speedup. Actual run receipts,
rather than the implemented workflow, determine which groups were measured.

The next exact-fixture adapter is `isogeny-perf kernels`, documented in
`KERNEL_BASELINE_PROTOCOL.md`: affine/cyclic Vélu, two x-only variants, two
square-root variants, and supplied-kernel Kohel. Its 84-row matrix is separate
from the 96-row BMSS series. It verifies degree/codomain and deterministic
point-transport fingerprints against both Kohel and Vélu; it does not claim
the BMSS series' full coefficient comparison. Two Apple M4 Pro runs now pass
84/84 paired fingerprints with all four counters identical; raw data and the
report are in `kernel-baseline-20261015/` and `kernel-report-20261015/`. Other
physical CPU baselines and kernel speedups remain unestablished.

CPU plan: Apple M4 Pro locally; GitHub native Linux x86-64, Linux ARM64,
Apple ARM64 and Intel macOS runners. Runner labels select architecture/OS,
not a guaranteed processor model. Each result records the actual CPU and
binary. [GitHub runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
