# Evidence status

Before outcomes: hash-only preregistration. The merged #781 base CNF is
archived once at `../rotated_s3_o_branch_20260925/evidence/producer/n13-m5/base.cnf`.
The frozen SHA-256 in `FROZEN.json` identifies it; no 96-copy base archive
is created. After host release, `smoke/` and `panel/` will retain raw
stdout/stderr, per-child JSON, exact derived-input hashes and commands.
CI independently reconstructs each derived input from the pinned base and
positive unit assumption, then replays output classifications and exact
signed point witnesses without starting any solver.
