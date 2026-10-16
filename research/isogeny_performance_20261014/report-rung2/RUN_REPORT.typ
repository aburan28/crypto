#set page(paper: "a4", margin: 18mm)
#set text(font: "Arial", size: 10pt)
= Native BMSS run record
2026-10-14 · Source-linked baseline and candidate

Requested: baseline each native isogeny-finding method across multiple processor types, then iterate on speed. This run implements the first eight-method series; the larger inventory remains open.

#image("PIPELINE.svg", width: 100%)

The baseline was frozen before algorithm edits. Deterministic supersingular fixtures span 192, 224, 256 and 511 field bits, at prime degrees 7, 31 and 101. These fixtures have their full coefficients and hashes in the raw JSONL; they are not the named P192/P224 models. Counts and native constructor verification run separately. Inversion internals are opaque.

The committed native comparison accepts all 96 exact output fingerprints and matches each input digest. Independent Kohel kernel and rational-map agreement and four public point transports guard every case. Missing terminal receipts, changed inputs and missing algorithm rows fail closed.

This machine is an Apple M4 Pro, 14 logical cores, 48 GiB, macOS arm64, compiler rustc 1.93.1 (Homebrew), as recorded by Cargo's compiler fingerprint. Cargo was invoked through rustup 1.98 but resolved the Homebrew rustc on PATH. It was shared with unrelated native jobs. No constructor wall-time sample is credited. All preparation, verification, counting, compilation and interrupted attempts are outside a claimed constructor interval, and no full-process speedup is claimed.

Failed attempts include system-disk ENOSPC, an interrupted queued build, an invalid revision argument, a 192-bit fixture using too many limbs, and incomplete sparse-checkout test fixtures. Each remains recorded in ATTEMPTS.md and the raw logs. Corrected runs use fresh paths. Compiler outputs are on the second SSD; evidence remains in the source worktree.

Comparison SHA256: #text(size: 7pt)[da9c0dd2e21396ad011cdeb4e23bfeb452f150df859689f2e1fb2b9ae157884d]

See REPORT.pdf for arithmetic counts; PROTOCOL.md for requirements; INVENTORY.md for the remaining methods; and the JSONL/COMPARISON.json files for raw inputs, binary hashes and exact output matches. CPU matrix measurements are credited only when their actual receipts are available.
