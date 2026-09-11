# Stage 28: independent internal relation-yield replay

Stage 28 closes the three internal payload-replay items left open by Stage 21. A separate fixed-width implementation reconstructs the `F_(2^23)` field, the binary Koblitz curve, the degree-12 `[0,2]` linearized-kernel factor base, the ordered rational points, every canonical pair sum, and every retained point witness. It imports no Stage-21 curve or factor-base code.

The replay regenerated 4,096 factor-base abscissae and 4,281 rational points. Its ordered factor-base BLAKE3 is `9b5bb635f1c505cb871a9e8723d2b79109688606cc42323ab92e05b8beaf4b58`. It enumerated all 9,165,621 pairs with `i <= j`, obtained 5,575,848 unique target entries and 3,589,773 duplicate pair sums, and reproduced canonical transcript BLAKE3 `c7a8ab74fc4ed9ea3832244f99b92e87b75ea268ffbcc89fb72fc15cd069ccae`.

All 384 retained target rows agreed with the independently built table. Exact curve re-addition verified 227 witnesses: 163 of 256 natural targets and all 64 planted targets. All 64 exact-miss controls remained absent. The metered replay used 7.999221 core-seconds, 8.442785 wall-seconds, and 146,636,800 bytes peak RSS.

This completes the finite public-synthetic internal payload replay and supports admitting the narrow Stage-21 relation-yield measurement within the project. It is not unaffiliated external reproduction, a SAT performance result, an end-to-end index-calculus experiment, or a Koblitz SOTA result.
