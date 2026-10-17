# Preserve replay validity when a Tier 1 host is noisy

The frozen dispatch and producer run `39048976945` remain unchanged. Its
four archives all recorded `PASS_INDEPENDENT_REPLAY`; strata 0 and 1 returned
`BOUNDED_UNKNOWN` because calibration exceeded the 10% gate. The original
post-run aggregator treated those expected censored strata as
`VERIFICATION_FAILED` and required their provider jobs to conclude
`success`, although the cloud verifier deliberately exited 2 for a noisy
host. The analysis repair separates the independent replay result from
host admission. It accepts provider `failure` only when that stratum's
archived analysis is `BOUNDED_UNKNOWN` and its independent local verifier
also exits 2. The full selector still requires all four `PASS_TIER1`
strata and therefore returns `HOST_STRATA_GATE`, `BOUNDED_UNKNOWN` for this
run.

| File | SHA-256 frozen in dispatch | Repaired SHA-256 |
| --- | --- | --- |
| `aggregate.py` | `9f0face7a9f4e0c336477c195ba7961b5f872341fef9c10bec3339845ec8a0f0` | `78fc24e6c45490edcfc696e8c9cd54a8ec5446752671dbfc8a9aa407b4511238` |
| `provider_jobs.py` | `d3b52618f84dae71c4d091a86544585ceecba35ce23b46b73559fe7b104c5c4e` | `10199210d23dfef8f60b0fd67c04b807517338fa42492335e0de94cf14628403` |

A direct 24-point aggregate check over the versioned cloud analyses produced
`PASS_INDEPENDENT_REPLAY` at the evidence layer and
`HOST_STRATA_GATE`/`BOUNDED_UNKNOWN` at the selector. The original raw ZIPs,
provider results, and calibration readings remain immutable. Full archive
binding awaits local replay of the two noisy strata; the admitted strata
have independent local replay receipts.
