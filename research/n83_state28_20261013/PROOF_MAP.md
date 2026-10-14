# N83 state storage proof map

| Result | Lean declaration | Hypotheses | Status |
| --- | --- | --- | --- |
| An 83-bit root survives an 11-byte modulus | `N83State28.root83_fits_eleven_bytes` | `x < 2^83` | Built; standard `propext` only |
| Packed state-array bytes under the source count bound | `N83State28.k2202_packed_state_bytes_bounded` | `states <= 83*2202^2`, 28-byte element | Built; no axioms |
| Maximum-size 48-to-28-byte replacement saves over 7 GiB | `N83State28.k2202_maximum_saves_more_than_seven_gib` | Numeric maximum, 48/28-byte element sizes | Built; no axioms |

The Rust `size_of` test checks the actual 48/28-byte layouts on the build
platform. The Rust round-trip test covers boundary 83-bit roots and both
stored root positions. The complete source-loop state bound, Rust byte-array
encoding bridge and first-witness relation identity are not definitions in
this Lean module; source inspection and the frozen native replay address
those links. The stronger end-to-end mathematical formalization remains open.
The pinned Lean 4 v4.34.1 build and exact `#print axioms` output are in
`proof_check.json` and `proof_check.log`. `lean4checker --fresh` was not
available on this host.
