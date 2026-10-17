# Prime-degree isogeny world-frontier campaign

## Current checkpoint — 2026-10-16

The [dated comparison and campaign checkpoint](record-comparison-20261016/README.md)
links the terminal expanded-map attempt, the separate compact result, the
fixed-modulus implementation controls and the remaining record-comparison gaps.
The [dedicated executable is officially published](https://github.com/aburan28/crypto/releases/tag/v0.1.2390);
use the [installation guide](../../docs/challenges/isogeny-prime-degree/README.md#run-the-dedicated-executable).
Publication evidence is retained in
[the release QA report](replay-optimization-v1/official-release-v0.1.2390/README.md).

The sections below preserve the initial planning and execution checkpoint.
References there to a running replay, a proposed compact category or pending
publication describe that earlier time; they are not current run status.
No frozen map, binary, season or earlier observation has been replaced.

## Initial campaign checkpoint

This active goal pursues the largest independently verified prime-degree
isogenies, with every registered standard prime-field model ranked separately.
The main challenge requires an expanded rational map and monic kernel. A second
research category will investigate compact evaluation routes without mixing its
scores with expanded-map results. The Codex runtime goal is anchored to chat
`01a11ebd-9436-73b3-b24e-63b10aad80fa`; its frozen harness goal is GOAL-ISOU-f30a07. The protocol, ownership and native
observations are published in [Crypto Autoresearcher PR #2221](https://github.com/aburan28/crypto-autoresearcher/pull/2221).

## Requirements and present evidence

| Requested requirement | State | Evidence or next gate |
|---|---|---|
| Work toward a world's-largest prime-degree isogeny | Active, no record claimed | Primary-source survey below; category and comparability must be established |
| All registered standard prime curves ranked separately | Existing implementation verified | Merged PR #1940; 54 exact standard models, 47 adapters, seven visible conversion gaps |
| Exceed the P192 degree-10453 incumbent | Construction completed; independent replay running | P192/50341; RUN-ISOU-f9546d, not admitted yet |
| Compact evaluation research | Survey only | Separate representation contract and independent evaluator required before execution |
| Reproducible artifacts and complete costs | Required per attempt | Frozen binary, catalogue, source identity, seed, native receipts and all failed attempts |
| Durable reviewable work | In progress | This branch and a separate harness protocol PR |

## Primary-source survey, 2026-10-15 PDT

This is a lower-bound survey, not an exhaustive record register. None of these
external calculations was independently reproduced in this round.

- **Jao and Soukharev (2010), sections 5.1–5.2.** Section 5.1 evaluates a
  degree-5000000029 isogeny over the field with p=10000000019 and gives the
  codomain, public-point image and conjugate-isogeny consistency check. Section
  5.2 evaluates the least Elkies prime greater than p/2 on Certicom ECCp-109,
  where p=564538252084441556247016902735257, again giving a codomain and point
  image. These are ideal-based evaluations; the paper does not publish expanded
  degree-sized rational-map coefficients. The ECCp-239 example in section 5.3
  records an ideal relation and explicitly leaves evaluation unfinished.
  [Versioned paper](https://arxiv.org/pdf/1002.4228v2).
- **Bernstein, De Feo, Leroux and Smith (2020), introduction.** Square-root
  Vélu evaluates quotient isogenies from a supplied prime-order kernel generator
  using softly square-root-many field operations. The paper distinguishes this
  task from emitting all coefficients of the defining rational functions.
  [Paper](https://arxiv.org/pdf/2003.10118).
- **Parallel square-root Vélu (2024), section 4.** The authors benchmark
  construction and two point evaluations for degrees from 19 through 321193
  over a chosen 1792-bit prime field. These are chosen-curve evaluations, not
  expanded maps from our standard sources. [Author preprint](https://www.preprints.org/manuscript/202401.1366),
  [implementation](https://github.com/TheParallelSquarerootVeluTeam/Parallel-Squareroot-Velu).
- **PRISM (2025–2026).** The authors provide implementations using compact
  representations of large prime-degree isogenies. Their degree parameters,
  representation and verification need a dedicated comparison before this
  campaign can name a global evaluated-degree target. [Author paper](https://r98inver.github.io/papers/2501_prism.pdf),
  [replication dataset](https://doi.org/10.48804/FZSO6R).

The surveyed evaluations already exceed degree 50341. A successful first rung
would advance the repository's P192 expanded-map frontier; it would not establish
a worldwide evaluated-degree record. The precise largest comparable published
expanded map is still unresolved. Arbitrarily large mathematical degrees,
inseparable Frobenius maps, composite-chain products and unevaluated ideal
relations do not discharge the prime-degree construction goal.

## First bounded execution

The Coordinator-approved protocol and execution card were committed and pushed
before the native trial. The installed executable runs this command:

```sh
isogeny challenge run --curve p192 --ell 50341 --out FRESH_ATTEMPT
```

Controls: native P192/13 and P521/7 positives, nine adversarial entries including
updated-hash kernel and map mutations, and fresh separate-curve leaderboard replay. The high-degree attempt retains the built-in
1800-second constructor watchdog, 3600-second replay watchdog, 8 GiB sampled RSS
limit and 256 MiB map limit. Those limits protect the machine; the campaign has
no fixed batch ceiling. A timeout is an incomplete attempt, not a statement of
nonexistence. One verified eigenline is sufficient, with coverage explicit.
The original relative-path control RUN-ISOU-ea3515 failed at its final leaderboard
assertion; its raw failure is preserved at harness commit b05e3aeca8. An additive
approved absolute-path amendment was published at 2e3220d223. The fresh control
RUN-ISOU-13a6b5 passed in full, then RUN-ISOU-f9546d started on 2026-10-16 at
06:53:59 UTC. The constructor emitted a 7,656,861-byte artifact, SHA256
82e3951488c793c5901ecf1dfeb3ab072b7cda16730899330df7caa9832a170c,
and returned PASS. Independent replay is running; no degree is admitted from
these interim observations. The finished control is published at 6d9015f199.

All attempts write directly into durable run directories. Every terminal round
gets raw receipts, a source-linked report, a diagram and inspected PDF. New
canonical curve edges or leaderboard values require independent replay and the
appropriate review; existing evidence is never overwritten.

The IC scoreboard, progress timeline and performance charts carry no new values
from this planning/survey round: it establishes neither an IC/rho ratio nor a
speedup. The existing challenge board retains its original pinned season until
new submissions are independently admitted under a specified season.

### Representation check: PRISM variants

The 2025 PRISM paper's basic response has composite degree q(2^a-q), even
though q itself is prime. A large prime factor alone does not earn a prime-degree
score. The 2026 journal extension, section 5, separately defines a prime-degree
variant with q=2^(2e-1)-s^2 and a higher-dimensional compact representation.
Its reference to that variant is not an independently replayed largest numeric
q. Auditing its implementation, extracting a reproducible instance, certifying
q and independently checking its evaluation are outstanding work.
[Journal article, sections 3 and 5](https://link.springer.com/article/10.1007/s00145-026-09585-6).
The [published C code](https://github.com/mariascrs/PRISM_v2) identifies itself as
the salted variant; it must not automatically be treated as the section-5
prime-degree variant. None of this code was executed in this survey round.

## Frozen native inputs

The first batch uses the tested CI merge executable at source
`6547d252889ac0e1089a0497bf80f05b884ab24e`, SHA256
`8e0ab524008d41098a816d117d563bbefa456959bb790a78e41efc41b05e6cc1`.
Its immutable archive is
`store:evidence/isogeny-world-frontier-20261016/isogeny-6547d252-macos-arm64`.
The existing compiled native admission controls are separately archived as
`store:evidence/isogeny-world-frontier-20261016/native-challenge-controls-macos-arm64`,
SHA256 `ccf19b8fdb15fb055059f9e2db7412a11d44fa526b080ae09462dc659debb527`.
Their custody receipts and the tracked `store/index.json` bind bytes and locations.

The frozen policy/catalogue season is
`f045887ec81e41c358ea3fa6f6bf15cfb33b8ed199e6ff589ffad66e8c60b457`.
See `frozen-policy.json` and `frozen-version.json`. All scientific trials use this
same executable even if main advances. The official ReleaseMe workflow 38857826201 completed every build job, including
all four isogeny platforms, but its first publication attempt failed while
pushing tag v0.1.2356 with HTTP 500. A failed-job-only retry was requested on
2026-10-16 at approximately 06:50 UTC. Publication remains unverified and its
binaries are not substituted into this frozen batch.

The initial report is one page and was rendered with Typst 0.15.1 and inspected
as a full-page image. The host has an Apple M4 Pro, 14 logical CPUs and 48 GiB RAM,
running Darwin 25.6.0 arm64. It is shared, with no isolated speed claim. The source-base release-library gate completed: 4148 passed, zero failed,
117 ignored; terminal exit code 0. It ran from the unchanged root source at
05e414f3fd908bde38f0e233b4d850270c025b07 with locked offline dependencies,
release optimization and one test thread. Root source and Cargo files are
unchanged by this study. Full terminal output and exit status are retained in
root-release-retry2.log and root-release-retry2.exit-code.txt. Earlier
interrupted validation attempts remain in validation-interruption.md and their
raw logs. The finished installed native control suite also passed in full at
harness run RUN-ISOU-13a6b5. No examples or Python implementations changed.
