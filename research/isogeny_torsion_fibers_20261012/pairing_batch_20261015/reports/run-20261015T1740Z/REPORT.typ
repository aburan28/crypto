#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= V6 batch pairing versus direct prime-field rho
*Run:* 2026-10-15 17:40 UTC. *Host:* shared Darwin arm64 L0, with concurrent release tests. *Spec:* `RUN_SPEC_V6.json`.

On 66 distinct order-67 targets, each of three arms returned the same verified scalars. Independent source-curve arithmetic checked 594/594 outputs. Complete cold CPU includes public-input charge, pairing-field construction and H search, reusable finite-field BSGS, every target solve and check, plus rho restarts. Three rotated rounds used the same inputs and caps.

At N=64, median complete CPU was 0.650 ms for direct F401 rho, 1.413 ms for extension-field rho, and 1.547 ms for pairing. Paired pairing/direct-rho ratios were 2.523, 1.620 and 3.000. The direct rho control was faster in all three rounds. Extension rho was a weaker implementation control because the source subgroup is F401 rational.

#image("GRAPH.svg", width: 100%)

These are shared-host L0 results. No counted-unit S or isolated scaling conclusion is available. Canonical frontier unchanged. The raw combined records SHA-256 is `fef2ab3733bd1912c2de0a144a8dd83f4f7c6563544db4ad6a78e70ed1742f24`; the independent replay receipt SHA-256 is `af26814d85377ef24487818b61571cd725972f4563bbf312f4fee825dc0775e9`. See `REPORT.md` for the exact protocol, links, limits, and older version receipts.
