# Stage 193 protocol: current native WDSat and CryptoMiniSat comparators

## Purpose and claim boundary

Run the two locally available named SAT comparators on the exact public
`n=59, ell=9, m=3` instance used by selected F4 and direct MITM. Replace the
historical externally-killed CryptoMiniSat receipt with an internally limited
run that can emit final conflict statistics, and independently rebuild the
exact configured WDSat source instead of trusting a stale binary.

This is a same-instance decomposition comparator. A timeout or solver UNKNOWN
is censored and never called `UNSAT`. It does not complete licensed Magma,
relation yield, an unknown-scalar DLP, full rho, external reproduction, novelty,
or SOTA.

## Frozen inputs and tools

- Crypto repository base:
  `637bc4d7cf3b600517c5bca823c8d8d026613dfa`.
- Source manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, with target-subgroup
  enumeration false and discrete-log-label use false.
- WDSat ANF SHA-256:
  `34c516d2568e327a3c5b403d6871845c80d523f0b7f012bd183be0bbc1469508`.
- CryptoMiniSat XOR-DIMACS SHA-256:
  `a529c4d8bb1e1a4b5904134b1e5d099adfb80806215f1ff128088d7d518d4bb8`.
- WDSat repository: `https://github.com/mtrimoska/WDSat` from the existing
  `/Volumes/SSD990/WDSat` checkout.
- WDSat source commit:
  `55d55b2620d768d9f7c78dcd8990a0689533c1d0`. This is the commit recorded by
  the sealed compatible build; `62be8cf6...` is the historical crypto
  repository Phase-B implementation commit and is not a WDSat object.
- WDSat configured `config.h` SHA-256:
  `4a73c3b5a14ded98f749d282b597594ae3ac05355a39cadcb9345a90bf735352`.
- WDSat static capacities: ANF variables 79, degree 4, internal variables
  2000, buffer 32804, OR equations 6511, XOR equations 112, XOR width 2001.
- CryptoMiniSat executable: `/opt/homebrew/bin/cryptominisat5`, version
  `5.14.7`, observed executable SHA-256
  `a3f85c3709b5e2a040bf82a4a604d1c7b9f10219bbf180a9e0f72319a2e892ac`.
- Native process meter: the Stage 192 selected binary at
  `/Volumes/SSD990/kic-stage192-target-e575d53f1/release/examples/koblitz_f4_stage188`,
  SHA-256
  `971ddd3b0a77b86aa1388632a201747cec42d570e4d9ba9536c7c61e8361c524`.

Create the detached WDSat worktree
`/Volumes/SSD990/wdsat-stage193-55d55b2` at the frozen source commit, install
the exact configured header, and run separate clean and native-build commands
under the Rust process meter. Rehash all 16 non-generated source and header
files against the sealed compatible source inventory. Record all
source/config/binary/input hashes. The repository-root WDSat binary and its
incompatible 52-variable configuration are forbidden.

## Frozen executions

All external thread caps are one. Use fresh process groups and native `wait4`
accounting.

### WDSat

Copy the exact ANF bytes to a short task-specific path to avoid the known WDSat
fixed-path-buffer failure. Run the freshly rebuilt binary with:

```text
-i /tmp/.../instance.anf
-g 1,2,...,27
```

Use a 120-second outer watchdog. A complete `SAT` or `UNSAT` terminal may be
accepted only with its printed conflict count. A watchdog kill without a
terminal is `timeout_inconclusive`, conflicts `null`.

### CryptoMiniSat

Run version 5.14.7 on the exact XOR-DIMACS input with:

```text
--threads 1 --verb 1 --verbstat 1 --printsol 0
--zero-exit-status --maxtime 120
```

Use a 135-second outer watchdog. Parse the final `s` terminal and `c conflicts`
statistic. `INDETERMINATE` is `timeout_inconclusive`, with the printed conflict
count retained. A SAT model would still require independent source-system and
curve-group validation before acceptance.

## Verification and accounting

A Rust-native Stage 193 verifier must:

1. authenticate input, tool, source, configuration and binary identities;
2. replay process-receipt hashes, return codes, watchdogs, wall, total CPU and
   peak RSS;
3. parse terminal state and conflicts without interpreting missing output as a
   proof;
4. bind both inputs to the same manifest/source instance and factor-base
   contract; and
5. include the verified Stage 192 F4/direct boundary by exact result hash only,
   without double-counting it as new Stage 193 work.

Charge WDSat worktree/config preparation, clean build, exact tests if any, both
solver processes, failures, verifier build/execution, composition, all thread
resources, memory and wall time. Complete campaign cost remains `null` unless
all inherited and new work is measured.
