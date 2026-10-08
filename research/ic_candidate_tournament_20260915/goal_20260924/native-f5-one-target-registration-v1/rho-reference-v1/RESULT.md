# Same-point strong-rho observation: exploratory only

The one declared `StrongRho` run recovered scalar `14668` for the exact public
point `[61889,74818]` used by the audited F5 target control. Its polynomial
basis scalar replay passed inside the timed interval, and the independent
checked-Sage certificate in `../result-v1/sage-replay.json` confirms the same
`14668·G = Q` equation. The original compact result and empty stderr are
retained as `rho-original.json` and `rho-original.stderr`.

| Record | Value |
| --- | ---: |
| Target validation and normal-basis conversion | 4,000 ns |
| Strong-rho walk and collision recovery | 36,625 ns |
| Independent polynomial-basis scalar replay | 15,334 ns |
| **One-target online interval** | **55,959 ns** |
| Target-independent curve and jump setup | 433,875 ns, separate |
| Walk steps / walks started | 88 / 40 |
| Distinguished-point entries at stop | 8; peak RSS unknown |

The three online phases sum exactly to 55,959 ns. The walk used signed
Frobenius, 32 lockstep lanes, four distinguished-point bits, one target and
one worker. The jump/start seeds were 2026100502/2026100503. The raw
operation counters in `rho-original.json` include 32 target-independent jump
scalar multiplications; no normalized operation score or memory peak is
inferred from them. No cross-target table or prior target collision work was
used. The source returned `verified_recovery=true` with no output on stderr.

The release example's source SHA-256 is
`b7aa13e342d87a7f17cf5be8b609cb665e1daaf8a56f7572a5eb7597e93155fc`;
its Cargo.lock SHA-256 is
`f99127c279e83c5fd8a474a604e5457a96d0ed16fb4c84758c78f51f1291ce79`.
The binary SHA-256 was
`93a07f4b282d30f5087fb9494e1bf82abad2d7f25fe44f43d894975896d8d09d`
both before and after the run. The example source first entered committed
revision `818e522ae`; the subsequent protocol-only revision `40c3f6a24`
changed no `src`, `examples`, Cargo or build source files. Rust/Cargo 1.93.1
built the example offline in release mode on macOS ARM64. The original result
SHA-256 is
`02e9987b60a8bbc27b8fd37b5364c6e5897858fbc5fbf063d722170088f47adf`.

For context, the separately frozen F5 worker's source-attested online interval
on this same point was 6,501.300958 ms. These are **diagnostic observations**:
the Mac had no auditable host isolation/noise receipt, the arms were not
interleaved, the rho binary lacks its own one-use frozen controller, the point
is not globally certified fresh, and the campaign's qualified
`rho_pairinv_4` has not been run on it. The row correctly retains
`source_bound_execution_admitted=false`, `fresh_paired_qualification=false`,
`headline_eligible=false`, null canonical IDs and `online_speedup=null`.
Neither an arithmetic ratio of these wall observations nor a global IC/rho
winner is claimed.
