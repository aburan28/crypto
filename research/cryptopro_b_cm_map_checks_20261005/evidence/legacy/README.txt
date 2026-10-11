CryptoPro-B public CM map correctness checks

Run python check_cm.py. Python standard library only; no third-party dependency.
The script evaluates the supplied degree-5 and degree-31 maps and final
isomorphism on known generated points of the actual CryptoPro-B curve.

Checks: component polynomial curve identities, source/target chain, kernel
square denominators, final isomorphism, generator order, CM discriminant
identity, subgroup scalar action, omega^2-omega+[155]=0, additivity, negation,
and the proposed polarization matrix determinant on the rational subgroup.

All point inputs are generated from known test scalars. There is no unknown
scalar recovery, relation collection, cryptanalytic parameter sweep, or solver.
The independent affine control checks forward-map agreement on these inputs.

This verifies the elliptic CM endomorphism. It does not verify an explicit
genus-two transfer: the supplied construction has no genus-two equation or
transfer formulas. It is not a relation-solving benchmark or a speedup claim.

The map file comes unchanged from the current cryptopro_b_evidence.zip;
its hash is recorded in results.json for reproduction.
