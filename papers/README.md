# Papers

Every manuscript generated from this repository's results lives here, one
directory per paper.

| directory | title | status |
|---|---|---|
| [`structural-completeness/`](structural-completeness/) | Structural Completeness of Prime-Field Elliptic Curves Against Isogeny-Graph Cryptanalysis | research draft, anonymized, not submitted; Markdown twin in `research/notes/cm-isogeny/PAPER_STRUCTURAL_COMPLETENESS.md` |
| [`lll-degeneracy-note/`](lll-degeneracy-note/) | An Empirical Note on LLL Degeneracy of the Boneh–Venkatesan Lattice for secp256k1 | research note, anonymized, not submitted |

Markdown drafts that are not yet typeset stay with their notes and are listed
here so they are not lost: `research/ecdlp_autolab/PAPER_DRAFT.md` (An Empirical
Resistance Map for Prime-Field ECDLP) and `research/iso1_weak_classes_20261007/REPORT.tex`
(a report, not a paper). Move one here when it becomes a manuscript.

Conventions shared by every paper here and by the sibling repositories
(`aburan28/cryptanalysis`, `aburan28/crypto`, `aburan28/crypto-autoresearcher`,
`aburan28/ml-cryptanalysis`), as set out in the "Papers" section of `AGENTS.md`:

- one directory per paper, holding the root `.tex` (or the `.md` of a draft not
  yet typeset), its sections and bibliography, a `README.md` (status banner,
  what it says, what it does not say, how to build, reproduction pointers, the
  checklist owed before submission), an `abstract.txt` and the built PDF;
- the build is `paper.mk` in this directory: a paper's `Makefile` sets `TEX`,
  `SECTIONS` and `BIB` and includes it (`make`, `make check`, `make arxiv`,
  `make eprint`, `make deanonymize`, `make clean`);
- every number cites a committed artefact in this repository or a sibling
  (`crypto:`, `cryptanalysis:` or `crypto-autoresearcher:` prefix); a stage
  diagnostic is never reported as a speedup; every end-to-end row has verified
  answers; every change is classed as advance, engineering, relabelling or
  accounting;
- `make check` must pass (references and citations resolved, no LaTeX errors,
  no overfull boxes) before a PDF is recommitted, and the PDF is rebuilt
  whenever a cited number, conclusion or figure changes;
- sources are anonymized (`\author{Anonymous}`); `make deanonymize` restores
  the authors from an untracked `AUTHORS` file copied from `AUTHORS.example`.

Add a new paper as a new subdirectory and a row in the table above.
