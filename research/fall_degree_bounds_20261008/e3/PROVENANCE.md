# E3 binary provenance

- `run.sh` built `hky_last_fall_degree` at commit `80e59d18`. Binary sha256
  `6c6a4a72071c2fe081520a3703e41ef70be1d4428c8b1675cfabe9b576892b7c` (`sha256`).
  The cells at `n = 18, 20, 22, 24, 26` ran on it.
- At 2026-10-12T00:36:33Z, while E3 was running, the binary was rebuilt in place to add
  the E4 `ladder` mode and an early stop in the closure ("`1 ∈ V` ⇒ `V_{F,c}` = all of
  `R_{≤c}`"). The rebuilt binary's sha256 is
  `fddf784d87302a8481456beb2d10456373d86b6fc4542cf3421c2301e9e8b6ac`.
- The cells at `n = 28` (started 00:36:39Z and 00:36:42Z) and `n = 32` ran on the rebuilt
  binary.
- **Why outputs are the same.** The early stop fires only once `1 ∈ V_{F,c}`. At that point
  the closure is the whole space, so the reported dimension is the one the full closure
  would reach. The `replay` check against the 32 Sage subfield records was rerun on the
  rebuilt binary, and all agree (`../e5/replay_subfield_rebuilt.jsonl`).
