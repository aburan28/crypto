# Manuscript validation — 2026-10-17

- Final pdflatex/BibTeX build: 12 pages. `make check` passes:
  all references resolved, no undefined citations, no LaTeX errors,
  no overfull boxes.
- All final PDF pages rendered and inspected. The two revised appendix pages
  were re-rendered and inspected after the citation correction.
- Both editable figure sources and PDF/SVG/150-dpi PNG outputs reviewed.
  Benchmark PDF and SVG/PNG copies are byte-identical to canonical files.
- Fifteen existing public evidence files were read at immutable commit
  4619bb97f63d960d40ffce127154b4f54195135d and match their local bytes.
  The full E0 input is unpublished and identified by exact SHA-256.
  The E1e support note was unpublished at that commit and is copied unchanged
  into `data/support_predicate_proof.md` with its original source/hash retained.
- Every local evidence pointer was SHA-256 checked, including the full input.
  The four copied raw JSON records match their frozen originals exactly.
- Self-review checked coefficient expansion, constant term square-root step,
  nonzero-z sufficiency, even-degree AS handling, empty-domain/constant-H
  cases, repeated-source orbits, exact count nullability, and odd-field scope.
  This is not an independent mathematical or human review.
- No new solver experiment, performance benchmark, root-wide test suite,
  native shared library build, or Lean build was run.
- `conductor check --summary ... --scope path:papers/norm-rootfinding` failed
  because conductor is absent. No task reservation is claimed.
  `command -v lean` and `command -v lake` found no executable.
- Novelty clearance, Lean proofs, characteristic-aware identity certificates,
  original-system coverage, matched baselines, and the qualified speed gate
  remain open, as recorded in the paper and proof map.
