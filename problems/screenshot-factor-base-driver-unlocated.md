# Screenshot factor-base candidate driver is not located

- **Status:** open
- **Kind:** source identification and implementation gap
- **Affected work:** [prime-field baseline](../research/prime61_factor_base_20261016/README.md)
- **Next action owner:** current research lane; source pointer requested from the user
- **Last updated:** 2026-10-16

## Problem

The supplied screenshot describes an F83 problem with subgroup order 41,
`G=(26,58)`, `Q=(61,82)`, a supplied degree-7 isogeny, three inverse pairs,
and 1,000/2,000 candidates ranked over twelve paired runs. Its exact driver,
candidate distribution, map and operation-accounting definition have not
been located in the inspected default branch of `aburan28/crypto`.
The requested larger field is either GF(2^83) or p=2^61-1. A replacement
baseline at the latter modulus does not reproduce that candidate experiment.

## Evidence and uncertainty

- Observed: default-branch source pin
  `de0c9391480155ebdd5a8495d2c79e4508d4795d`; native `prime_fast` baseline
  and [frozen protocol](../research/prime61_factor_base_20261016/protocol.json).
  All twelve 61-bit runs exhausted their caps with zero relations; the
  20-bit positive control solved and verified.
- Observed: local source and repository code searches did not locate the
  screenshot's exact point pair, numeric averages or inverse-pair driver.
  Other branches, archives and earlier uncommitted sessions have not been
  exhaustively inspected. No claim of repository-wide absence is made.
- Derived or proved: no new mathematical obstruction established.
- Missing capability: original executable experiment contract and map.
- This evidence does not establish nonexistence of relations or a large-field
  isogeny, algorithmic infeasibility, or a candidate cost ranking.

## What would unblock it

- [ ] Identify the screenshot driver's file path, branch, PR or source bundle;
  recover its original target, map, candidate policy, seeds and counters.
- [ ] Preserve the algorithm in a native implementation, validate a fresh
  requested-field curve/subgroup/map, and freeze its larger-field run protocol.
- [ ] Run both original candidate-pool sizes with complete capped records,
  verified scalar recovery for any completed solve, and priced operation counts.

The first condition permits a faithful port; it is not itself a completed
large-field benchmark. A mathematical map and subgroup must be checked again
after a field change. Any binary-field version also needs binary arithmetic.

## Next discriminating step

Locate the original source pointer and replay its F83 fixture unchanged as a
control. Then preregister the requested-field port with the 1,000/2,000 pools
and twelve paired seeds. Record required resources from that concrete code;
no new large campaign is scheduled by this problem record. A missing source
pointer leaves this requirement open. A control mismatch must be corrected
before interpreting larger-field costs.

## Resolution and history

2026-10-16: bounded native prime-field baseline completed; screenshot
candidate search remains unresolved. Its report and this record do not
discharge the original request.
