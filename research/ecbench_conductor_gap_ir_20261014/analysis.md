# Complete implementation instruction census

Unit: Callgrind simulated user-space instructions, inside `methods::solve`. External independent verification is outside this interval. No native wall-time claim.

Verified profiles: 100/100. A/A maximum relative deviation: 0.002238%. Frozen repeated instruction-lead criterion: false.

| Curve | Arm | n | Mean Ir | Mean Ir/sqrt(r) | Min Ir | Max Ir |
|---|---|---:|---:|---:|---:|---:|
| icv1-fp29-t47-14392892 | rho-neg | 8 | 14656219 | 819.189790 | 4558199 | 29817057 |
| icv1-fp29-t47-14392892 | rho-neg-aa | 8 | 14656219 | 819.189790 | 4558301 | 29817057 |
| icv1-fp29-t47-14392892 | ic-mitm | 8 | 412609741 | 23062.271080 | 219310496 | 639986192 |
| icv1-fp29-t47-1f934956 | rho-neg | 8 | 15814761 | 883.944964 | 5582155 | 28297388 |
| icv1-fp29-t47-1f934956 | rho-neg-aa | 8 | 15814748 | 883.944251 | 5582053 | 28297388 |
| icv1-fp29-t47-1f934956 | ic-mitm | 8 | 399739050 | 22342.880961 | 141671255 | 552318765 |
| icv1-fp37-t277-91fd6e29 | rho-neg | 8 | 274100852 | 959.306807 | 92812878 | 477452242 |
| icv1-fp37-t277-91fd6e29 | rho-neg-aa | 8 | 274100885 | 959.306924 | 92812878 | 477452242 |
| icv1-fp37-t277-91fd6e29 | ic-mitm | 8 | 13094204544 | 45827.510043 | 9526898259 | 16142432711 |
| icv1-fp37-t277-090dc414 | rho-neg | 8 | 323334289 | 1131.615542 | 52948331 | 701693062 |
| icv1-fp37-t277-090dc414 | rho-neg-aa | 8 | 323334256 | 1131.615426 | 52948331 | 701693062 |
| icv1-fp37-t277-090dc414 | ic-mitm | 8 | 15222027071 | 53274.530432 | 8612452164 | 17440924750 |

| Comparison | Curve/gap | Ratio | 95% interval | 99% interval |
|---|---|---:|---|---|
| ic_over_rho | icv1-fp29-t47-14392892 | 28.152537 | 16.318500 to 56.635750 | 14.012291 to 74.394437 |
| ic_over_rho | icv1-fp29-t47-1f934956 | 25.276326 | 15.892807 to 40.760028 | 13.301184 to 48.247036 |
| floor_over_crater_ic | ell8209-q29 | 0.968807 | 0.657105 to 1.353537 | 0.573521 to 1.492155 |
| ic_over_rho | icv1-fp37-t277-91fd6e29 | 47.771484 | 36.856235 to 68.608686 | 34.841976 to 79.650822 |
| ic_over_rho | icv1-fp37-t277-090dc414 | 47.078295 | 32.205867 to 83.655392 | 29.449404 to 104.838149 |
| floor_over_crater_ic | ell131101-q37 | 1.162501 | 0.960003 to 1.385142 | 0.901657 to 1.462762 |

All individual profiles, source run/record IDs and validation failures are in analysis.json. Missing profiles invalidate inference; no successful-subset analysis is used. Intervals containing one do not prove equivalence.
