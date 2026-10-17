#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= Compact batch V5 versus extension-curve rho
*Run:* 2026-10-15 17:29 UTC. *Host:* shared Darwin arm64, L0, with concurrent release tests. *Spec:* `RUN_SPEC_V5.json`.

The compact native runner timed complete cold pairing and same-target extension-curve rho without per-target JSON output, adding public-input validation to both. It charged field construction, failed H search, alpha, reusable finite-field BSGS, each beta and solve, scalar checks and rho restarts. Separate arithmetic replay checked 108, 360 and 396 outputs on the order-19, order-61 and order-67 cells.

At prefixes `(p,r,N)=(61,19,16)`, `(271,61,60)` and `(401,67,64)`, median pairing / extension-rho CPU was 128.772 / 0.474 ms, 40.138 / 1.216 ms and 1.275 / 1.703 ms. Median paired ratios were 261.81, 33.60 and 0.889, respectively.

#image("GRAPH.svg", width: 100%)

The order-67 ratio fell below one in all three V5 rounds at N=64, but this rho arm computed over the extension field although the source points are prime-field rational. V6 adds that stronger prime-field reference. This V5 result is a representation-specific, shared-host observation, with no counted-unit or scaling claim. Raw order-67 SHA-256: `5c6a927cf5e1d6b570b80c009367e8cd98f1a84778c6b1ccd1441e6c8125d8eb`. See `REPORT.md` for all raw and replay links.
