#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
= Cached quintic Frobenius: theory cycle, 2026-10-17
*Status:* written derivation and prospective experiment, no executed benchmark.
The managed executor is offline. This source and the editable SVG are the
complete rendering input; no rendered or inspected PDF is claimed.

== Exact mechanism
In characteristic two, cache $X^6 mod H$ and $X^8 mod H$ for a monic quintic.
For $r=sum_(i=0)^4 r_i X^i$, square each coefficient and combine the cached
columns. This is Frobenius-semilinear over the field, and remains exact in a
nonreduced quotient. The GF4 example $omega^2=omega+1$ rules out a linear cache
that leaves coefficients unchanged. Odd characteristic needs a different map.

#figure(image("flow.svg", width: 100%), caption: [Proposed operator; no measured gain.])

== Conditional work model
With coefficient-square cost S and multiply cost M, the stated dense model
changes $105S+510M$ to $105S+315M$ for18 recurrence and3 tail steps, including
setup. Its ratio is at most34/21 under that model. This is not a wall-time bound
or an8x result. CONTRACT.md supplies all assumptions and the written derivation.

== Primary-source review
The repository's rootfind-mitm source stops at its first rational group-law hit.
It is not an all-original-polynomial-roots baseline. Generic wide-prime F4
backends now exist: do not describe all wide F4 arithmetic as missing.
Full matching extraction, exact field/domain identity and the fastest external
references remain separate obligations.

== Next experiment and remaining gates
The blocked backlog card XP-20261017-quintic-support-frobenius defines three
bounded native experiments, one CPU,4GiB and90seconds per worker. No task was
scheduled. Before readiness, pin code and input hashes, check independent
oracles, and validate the card. Require complete setup/remainder/extraction/
verification/serialization accounting, cold pairs and A/A for speed claims.
Preserve performance/GATE.md and the requested8x target.
Lean proofs, identity-certificate bridge, native replay, PDF rendering and
the complete matched hard comparison remain open.
