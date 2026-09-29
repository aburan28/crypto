| reference rho | Ir (whole process) | IC / rho | steps | verified |
|:--|--:|--:|--:|--:|
| R0 as merged in PR #955 (its own binary) | 371,102,176,689 | 0.240 | 19,103,507 | 1,024/1,024 |
| R0 regression: this binary, rho-local software field | 406,208,773,598 | 0.220 | 19,103,507 | 1,024/1,024 |
| R1 library Gf2 field (hardware clmul, table reduction) | 281,027,207,019 | 0.317 | 19,103,507 | 1,024/1,024 |
| R2 + normal-coordinate canonicalization, fast mulmod | 68,926,289,783 | 1.294 | 19,247,253 | 1,024/1,024 |
| R3 + 32 lockstep lanes, Gf2::batch_inv | 46,384,789,571 | 1.923 | 19,569,835 | 1,024/1,024 |
| IC `koblitz_orbit_dlp_fast` (current source) | 89,190,346,806 | — | — | 1,024/1,024 |

rho_best = R3 (46,384,789,571 Ir); rho* = IC/rho_best = 1.923 -> lead dies
