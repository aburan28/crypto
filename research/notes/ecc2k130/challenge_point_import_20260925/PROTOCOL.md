# Direct import of the public ECC2K-130 points into the rotated n=131 model

The n=131 support gate in [PR #762](https://github.com/aburan28/crypto/pull/762)
uses the bit-polynomial `z^131+z^13+z^2+z+1`. The established ECC2K-130
challenge implementation, [`relations.py`](https://github.com/aburan28/crypto/blob/a1df88fbf76fd14e565a5ec6b7d0ea4402f954ed/research/ecc2k130_relations/relations.py)
(Git blob `a1df88fbf76fd14e565a5ec6b7d0ea4402f954ed`), and its
[`fastfield.py`](https://github.com/aburan28/crypto/blob/bb57818dd85db7bff049e6f590d818616c5164e6/research/ecc2k130_relations/fastfield.py)
(Git blob `bb57818dd85db7bff049e6f590d818616c5164e6`) use **the same
polynomial and bit-i-is-z^i encoding**. Consequently the public P and Q
integers can be copied directly; a field-isomorphism computation is not
needed for these two repository models. The source is the public Certicom
challenge point pair, not a recovered scalar or a test target chosen for
decomposability.

The generator freezes the literal P/Q coordinates from `relations.py`, q,
cofactor four, and the lambda from PR #762. It requires canonical 131-bit
coordinates, `#E(F_(2^131))=4q` from the Weil trace recurrence, on-curve
P/Q, `[q]P=[q]Q=O`, `tau(P)=[lambda]P`, `tau(Q)=[lambda]Q`, and exact
rank-131 normal conjugates of beta 3. It emits polynomial and normal-basis
coordinates for P, Q, and the four exact `Q+T` torsion translates, with the
shared `[4]Q` projection. It checks the torsion kernel and covariance
`tau(Q+T)=tau(Q)+T`; multiplying a full translated point by lambda is **not**
the asserted identity.

An independent bit-serial/Fermat field and curve implementation rechecks
every fixture field coordinate, normal-basis mask, group identity, torsion
translate and public trace. CI regenerates the committed fixture byte for
byte and checks it with the independent implementation. The fixture is a
reusable literal input for future actual-point PDP exporter gates. It
contains no decomposition witness, relation rank, discrete logarithm,
runtime comparison, or n=131 support estimate. An actual-target solver
still needs complete direct/chain S6/S7 constraints, rational-lift and
exceptional-branch handling, a charged search and independent replay.
