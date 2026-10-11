# Point-decomposition preregistrations, 2026-10-07

Six protocols for the point-decomposition problem (PDP) in elliptic-curve
index calculus, frozen before any instrument is built or run.  Every row is
**PENDING**: nothing here is a measurement, a speedup, or a ledger change.
Each protocol states its derivation, its instrument, numbered pass/fail
predictions, its registered class, its stop condition and its inadmissible
moves, in the form of
[`../geometric_v_linear_20261006/PROTOCOL.md`](../geometric_v_linear_20261006/PROTOCOL.md)
and [`../graph_cycle_ic_20261007/PROTOCOL.md`](../graph_cycle_ic_20261007/PROTOCOL.md).
Instruments are Rust, per `AGENTS.md`; the Python instruments of the two
cited protocols are historical.

## The two filters every protocol was screened against

**Filter 1 — PDP is k-SUM on the curve group.**  Deciding whether a target
is a sum of `k` factor-base points using only group operations is k-SUM
over a black-box group: meet-in-the-middle costs `N^{⌈k/2⌉}` time, or
`N^{k−1}` time with `N` memory, and no generic algorithm is known to do
better.  The product law of
[`../notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md`](../notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)
§5 is this hardness in the harness's unit.  An idea has a chance only if it
uses the field representation in a way that is not enumeration of a
sub-tuple; every algebraic oracle measured so far (Gröbner, SAT, crossbred,
Riemann–Roch, Hamming ideals, the XL cores) was enumeration in disguise.

**Filter 2 — field additivity and group additivity align in one known
place.**  Subspaces are additive in the field; sumsets are additive in the
group; the `x`-map connects them only through `S₃`.  The single alignment
found is that squaring is `F₂`-linear, which makes `S₃` bilinear in the
bits of two summands and linear in a pair's symmetric functions
([`../geometric_v_linear_20261006/README.md`](../geometric_v_linear_20261006/README.md)).
D-1 and D-2 look for a second alignment or bound the first; D-3 tests the
one escape from the resulting gap formula; E-1 to E-3 move constants in the
compact-orbit producer without touching any exponent.

## The protocols

| id | protocol | what it would change | registered class | cost | kill line |
|:--|:--|:--|:--|:--|:--|
| D-1 | [non-even coordinates compatible with the Artin–Schreier form](PROTOCOL-D1-artin-schreier-coordinates.md) | the oracle's degree and sign fibre | boundary, engineering if it passes | days | per-block Boolean degree above the `x`-system's plus one |
| D-2 | [exact reach of linearised symmetric oracles at four summands](PROTOCOL-D2-linearisation-reach.md) | closes or opens the linearisation class with a number | boundary | 1–2 days | reach below 20 closes it |
| D-3 | [relation-matrix filtering against the linear-algebra cap](PROTOCOL-D3-filtering-la-cap.md) | whether a three-summand oracle could ever suffice | boundary | 1 week | filtering factor below 40 |
| E-1 | [packed root index](PROTOCOL-E1-packed-root-index.md) | bytes per state, hence `K` at fixed memory, hence probes per relation | engineering | days | at or above 64 bytes per state |
| E-2 | [batch verification](PROTOCOL-E2-batch-verification.md) | the largest charged non-algorithmic cost | engineering, accounting-preserving | days | any planted-bad column passes |
| E-3 | [rank from index self-collisions](PROTOCOL-E3-index-self-collisions.md) | replaces rank probing with a larger index | engineering or boundary | days | collisions off the law by more than 2× |

None of the six touches the `vs_rho` verdict of
[`../../docs/ic/BOUNDARY_TARGETS.md`](../../docs/ic/BOUNDARY_TARGETS.md): the
strong-rho ladder, the whole-process instruction counts and the frontier in
[`../../docs/bounds/FRONTIER.md`](../../docs/bounds/FRONTIER.md) stand until a
result lands under one of these protocols and is classed by §3 of
`AGENTS.md`.

## Requirement-to-evidence table

| requested | delivered in this PR | status |
|:--|:--|:--|
| D-1 to D-3 and E-1 to E-3 as preregistered protocols with kill tests | six `PROTOCOL-*.md` files, this index | drafted; every measurement PENDING |
| repo protocol format | derivation, instrument, predictions, stop condition, inadmissible moves, class, per protocol | done |
| one planning PR | this directory, no code, no ledger change | done |
| instruments | named per protocol, not built | not attempted here; each is a follow-on PR |
| Conductor task | `conductor check` was run and could not reach the local server (`Post https://127.0.0.1:8443/...: EOF`) | blocked; recorded, not bypassed |

## What this PR deliberately leaves out

These were recommended in the same review and are not preregistered here,
so that this PR stays one coherent thing:

- the preprocessing-class control for the online-after-precompute ladder
  (a Bernstein–Lange table on the curve at matched precompute, with
  `S·T²/N` and the Corrigan-Gibbs–Kogan bound as the floor);
- the counting note closing the ledger's "yield above ceiling" target for
  uniform targets;
- the isogeny-side protocols (class-level law for Joux–Vitse weakness,
  level-aware walks, per-curve rigidity certificates).

## Evidence cited from unmerged work

Three figures below come from branch `cursor/ic-boundary-experiments-d111`,
which is not on `main`: the `n = 83`, `a = 1` compact-orbit run
(`experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json`
at commit `7b1787074`: 29,880,000 regular states, 29,878,182 root-table
entries, peak RSS 5,771,837,440 bytes, `K = 600`, 8,845,441 target probes),
the zero-collision probe (`examples/koblitz_zero_collision_probe.rs`), and
the `u128` wide-path producer (`examples/koblitz_orbit_dlp_fast.rs`).  They
are labelled unmerged wherever they appear and are inputs to a protocol,
never evidence for a claim.

## Curves

Curves are named by ICV1 slug.  The toy and ladder curves these protocols
run on are already registered: `icv1-f2m13-t181-515ee569`,
`icv1-f2m19-t797-b6cf2467`, `icv1-f2m23-t5197-69e76b73`,
`icv1-f2m41-tm2308219-7f48b14a`, `icv1-f2m53-tm56619371-dac20a85`,
`icv1-f2m61-t158598901-ab42b6c5`, `icv1-f2m83-t6151469093347-cdcc5432`
and the §8a gate curve `icv1-f2m83-tm6151469093347-debefd74`.  The
challenge curve is ECC2K-130.  A protocol that needs a curve not in the
registry registers it in its own follow-on PR.
