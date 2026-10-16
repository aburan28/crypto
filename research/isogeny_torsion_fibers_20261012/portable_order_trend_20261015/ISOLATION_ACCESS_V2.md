# Isolated timing access follow-up, 2026-10-16

The [V1 access receipt](ISOLATION_ACCESS_V1.md) found no configured
`isolab`, a Standard On-Demand EC2 quota of zero vCPUs in five regions,
and an `AccessDeniedException` on the requested quota increase. This
read-only follow-up checked the account's **All Spot Instance Requests**
vCPU quota (`ec2`, `L-34B43A08`) in the same five regions:

```sh
aws service-quotas get-service-quota --region REGION \
  --service-code ec2 --quota-code L-34B43A08 \
  --query 'Quota.Value' --output text
```

| Region | Spot vCPU quota returned |
| --- | ---: |
| `us-west-2` | 0.0 |
| `us-east-1` | 0.0 |
| `us-east-2` | 0.0 |
| `us-west-1` | 0.0 |
| `eu-west-1` | 0.0 |

Each command exited zero. The screened 128-vCPU `c6i.metal` bare-metal
candidate cannot be placed through this account's Standard On-Demand
or Spot capacity under the returned quotas. No instance or quota
request was created in this follow-up. Other regions, accounts, host
types and authorized lab workers were not exhaustively surveyed, so
this is a placement result for the checked routes, not a universal
unavailability claim.

An account owner enabling sufficient capacity or an authorized
whole-core lab worker would reopen the preregistered interleaved L2
session. That session still needs CPU placement and interference
telemetry, a same-target rho A/A gate, full raw resource receipts and
independent scalar replay. Hosted GitHub VM CPU ratios remain L1
descriptions and do not meet those admission criteria. This access
limit is separate from the unresolved complete portable-price
instrumentation and from the mathematics of torsion fibers.
