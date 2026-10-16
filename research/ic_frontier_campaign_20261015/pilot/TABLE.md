# Phase 0 pilot: every record across the AGENTS.md §5a band

Produced by `ecbench band --ranks` (native) over `sessions/*`. Single train targets; not a measurement.

| curve | log2 r | arm | verified/runs | S(w=0) | S(w=0.1) | S(w=1) | uncounted |
|:--|--:|:--|--:|--:|--:|--:|:--|
| icv1-f2m13-t181-515ee569 | 10.97 | bsgs-neg | 2/2 | 1.531 | 1.668 | 2.905 | — |
| icv1-f2m13-t181-515ee569 | 10.97 | ic-compact | 2/2 | 14.747 | 24.640 | 113.675 | — |
| icv1-f2m13-t181-515ee569 | 10.97 | ic-f4 | 2/2 | 27.036 | 41810.222 | 417858.894 | — |
| icv1-f2m13-t181-515ee569 | 10.97 | ic-negfold | 0/2 | — | — | — | — |
| icv1-f2m13-t181-515ee569 | 10.97 | ic-shared | 2/2 | 56.586 | — | — | field_arithmetic_uncharged, hash_and_allocation_uncharged, modular_combination_uncharged |
| icv1-f2m13-t181-515ee569 | 10.97 | rho-frob | 2/2 | 1.743 | 1.773 | 2.044 | — |
| icv1-f2m13-t181-515ee569 | 10.97 | rho-strong | 2/2 | 13.972 | 14.135 | 15.603 | — |
| icv1-f2m19-t797-b6cf2467 | 17.00 | bsgs-neg | 1/1 | 1.039 | 1.140 | 2.048 | — |
| icv1-f2m19-t797-b6cf2467 | 17.00 | ic-compact | 1/1 | 3.118 | 5.233 | 24.264 | — |
| icv1-f2m19-t797-b6cf2467 | 17.00 | ic-f4 | 1/1 | 9.761 | 964887.235 | 9648784.502 | — |
| icv1-f2m19-t797-b6cf2467 | 17.00 | ic-negfold | 1/1 | 26.581 | 30.897 | 69.739 | — |
| icv1-f2m19-t797-b6cf2467 | 17.00 | ic-shared | 1/1 | 6.529 | — | — | field_arithmetic_uncharged, hash_and_allocation_uncharged, modular_combination_uncharged |
| icv1-f2m19-t797-b6cf2467 | 17.00 | rho-frob | 1/1 | 0.896 | 0.941 | 1.354 | — |
| icv1-f2m19-t797-b6cf2467 | 17.00 | rho-strong | 1/1 | 2.823 | 2.889 | 3.476 | — |
| icv1-fp20-t157-e39da206 | 19.07 | bsgs-neg | 2/2 | 0.584 | 0.641 | 1.152 | — |
| icv1-fp20-t157-e39da206 | 19.07 | ic-pair | 2/2 | 20.978 | 23.265 | 43.853 | — |
| icv1-fp20-t157-e39da206 | 19.07 | ic-rootfind | 2/2 | 1363.669 | 1500.786 | 2734.841 | — |
| icv1-fp20-t157-e39da206 | 19.07 | kangaroo | 2/2 | 2.790 | 2.819 | 3.085 | — |
| icv1-fp20-t157-e39da206 | 19.07 | rho-neg | 2/2 | 2.003 | 2.183 | 3.797 | — |
| icv1-f2m23-t5197-69e76b73 | 21.00 | bsgs-neg | 1/1 | 1.457 | 1.602 | 2.905 | — |
| icv1-f2m23-t5197-69e76b73 | 21.00 | ic-compact | 1/1 | 3.110 | 4.416 | 16.170 | — |
| icv1-f2m23-t5197-69e76b73 | 21.00 | ic-f4 | 1/1 | 16.599 | 7307108.058 | 73070931.181 | — |
| icv1-f2m23-t5197-69e76b73 | 21.00 | ic-negfold | 1/1 | 127.199 | 151.828 | 373.487 | — |
| icv1-f2m23-t5197-69e76b73 | 21.00 | ic-shared | 1/1 | 2.234 | — | — | field_arithmetic_uncharged, hash_and_allocation_uncharged, modular_combination_uncharged |
| icv1-f2m23-t5197-69e76b73 | 21.00 | rho-frob | 1/1 | 0.201 | 0.209 | 0.283 | — |
| icv1-f2m23-t5197-69e76b73 | 21.00 | rho-strong | 1/1 | 1.249 | 1.347 | 2.226 | — |
| icv1-fp24-t2705-27cbb01b | 23.38 | bsgs-neg | 1/1 | 0.522 | 0.574 | 1.040 | — |
| icv1-fp24-t2705-27cbb01b | 23.38 | ic-pair | 1/1 | 14.930 | 17.548 | 41.116 | — |
| icv1-fp24-t2705-27cbb01b | 23.38 | ic-pkm | 1/1 | 31.446 | 34.427 | 61.248 | — |
| icv1-fp24-t2705-27cbb01b | 23.38 | ic-rootfind | 1/1 | 1676.316 | 1845.199 | 3365.147 | — |
| icv1-fp24-t2705-27cbb01b | 23.38 | kangaroo | 1/1 | 3.974 | 4.000 | 4.232 | — |
| icv1-fp24-t2705-27cbb01b | 23.38 | rho-neg | 1/1 | 1.411 | 1.541 | 2.710 | — |
| icv1-fp28-t4155-50f492b1 | 26.20 | bsgs-neg | 1/1 | 0.515 | 0.566 | 1.028 | — |
| icv1-fp28-t4155-50f492b1 | 26.20 | ic-pair | 1/1 | 86.837 | 104.020 | 258.668 | — |
| icv1-fp28-t4155-50f492b1 | 26.20 | ic-pkm | 1/1 | 131.152 | 159.844 | 418.073 | — |
| icv1-fp28-t4155-50f492b1 | 26.20 | ic-rootfind | 1/1 | 11054.586 | 12168.592 | 22194.640 | — |
| icv1-fp28-t4155-50f492b1 | 26.20 | kangaroo | 1/1 | 2.686 | 2.703 | 2.861 | — |
| icv1-fp28-t4155-50f492b1 | 26.20 | rho-neg | 1/1 | 1.005 | 1.101 | 1.967 | — |
| icv1-fp32-t7319-90da43e9 | 32.00 | bsgs-neg | 1/1 | 1.106 | 1.217 | 2.212 | — |
| icv1-fp32-t7319-90da43e9 | 32.00 | ic-pair | 0/1 | — | — | — | — |
| icv1-fp32-t7319-90da43e9 | 32.00 | ic-rootfind | 0/1 | — | — | — | — |
| icv1-fp32-t7319-90da43e9 | 32.00 | kangaroo | 1/1 | 0.751 | 0.752 | 0.763 | — |
| icv1-fp32-t7319-90da43e9 | 32.00 | rho-neg | 1/1 | 0.987 | 1.085 | 1.968 | — |
| icv1-fp36-t454193-44a20818 | 36.00 | bsgs-neg | 1/1 | 1.413 | 1.554 | 2.826 | — |
| icv1-fp36-t454193-44a20818 | 36.00 | ic-pair | 0/1 | — | — | — | — |
| icv1-fp36-t454193-44a20818 | 36.00 | ic-rootfind | 0/1 | — | — | — | — |
| icv1-fp36-t454193-44a20818 | 36.00 | kangaroo | 1/1 | 3.016 | 3.018 | 3.039 | — |
| icv1-fp36-t454193-44a20818 | 36.00 | rho-neg | 1/1 | 1.344 | 1.478 | 2.687 | — |
| icv1-fp40-tm924137-0ecaaeb9 | 40.00 | bsgs-neg | 1/1 | 0.614 | 0.675 | 1.228 | — |
| icv1-fp40-tm924137-0ecaaeb9 | 40.00 | ic-pair | 0/1 | — | — | — | — |
| icv1-fp40-tm924137-0ecaaeb9 | 40.00 | ic-rootfind | 0/1 | — | — | — | — |
| icv1-fp40-tm924137-0ecaaeb9 | 40.00 | kangaroo | 1/1 | 3.382 | 3.384 | 3.396 | — |
| icv1-fp40-tm924137-0ecaaeb9 | 40.00 | rho-neg | 1/1 | 1.623 | 1.785 | 3.245 | — |
