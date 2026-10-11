# Isogeny-loop endomorphism sweep (Koshelev–Sanso relation-lattice GLV), 2026-10-10

Evidence for `RESEARCH_ISOGENY_LOOP_GLV.md`.  Every file here was produced by
`scripts/isogeny_loop_sweeper.py` at the revision of this commit, with cypari2
(PARI/GP 2.17.2) on the session container (4 cores, no GPU).  The `.md` files are the
tool's stdout, the `.json` files the full machine-readable result (all loops under
budget, dlogs, eigenvalues, GLV bases), `*.log` the stderr progress lines.

| file | command |
|:--|:--|
| `selftest_D-71.{md,json}` | `python3 scripts/isogeny_loop_sweeper.py selftest --json …` |
| `group_a_koblitz_and_known.{md,json}` | `python3 scripts/isogeny_loop_sweeper.py curve ECC2K-130 icv1-f2m83-tm6151469093347-debefd74 icv1-f2m83-t6151469093347-cdcc5432 ECC2K-95 ECC2K-108 secp256k1 id-GostR3410-2001-CryptoPro-B-ParamSet id-GostR3410-2001-CryptoPro-A-ParamSet sect113r1 --top 6 --json …` |
| `group_b_binary131.{md,json}` | `python3 scripts/isogeny_loop_sweeper.py curve ECC2-131 sect131r1 sect131r2 --top 6 --json …` |
| `group_c_certicom_mid.{md,json}` | `python3 scripts/isogeny_loop_sweeper.py curve ECCp-131 ECC2-109 ECCp-109 ECC2-97 ECCp-97 --top 6 --json …` |
| `group_d_certicom_small.{md,json}` | `python3 scripts/isogeny_loop_sweeper.py curve ECC2-89 ECCp-89 ECC2-79 ECCp-79 --top 6 --json …` |
| `scan_D3-3000_bits256.{md,json}` | `python3 scripts/isogeny_loop_sweeper.py scan --scan-from 3 --scan-to 3000 --bits 256 --top 4 --json …` |
| `scan_D3-3000_bits128.{md,json}` | same with `--bits 128` |
| `scan_D3-3000_bits128_char2.{md,json}` | same with `--bits 128 --char2` |

Scan and disc modes enumerate only up to the LLL bound (enough to certify the minimum) and list prime-norm loops only (the automorphisms of D = −3, −4 are not loops).

Defaults in force: primes ℓ ≤ 100, projective cost model (7.5ℓ M per ℓ-step, 2 M for
ℓ = 2 in characteristic 2), doubling 8 M, enumeration radius = the GLV saving 8⌈log₂r/2⌉,
`--max-disc-bits 140` (above it only the heuristic estimate is printed), maximal order.

Curve parameters: registry aliases resolve through `docs/curves/registry.json`; the
Certicom challenge curves that are not in the registry (ECCp-131/109/97/89/79,
ECC2-131/109/97/89/79, ECC2K-108) are built into the script from the Certicom ECC
challenge document (field, a, b, h, n in hex).

Conductor (the repository's task-coordination control plane) was not reachable in
the session that produced this directory; edits went through the hook's fallback.
