"""Exact postexecution negative proofs on bounded three-summand group inputs.

Every ordered triple is P_i+(P_j+P_k). The pair set includes repetitions,
opposite points and identity sums. Scanning Q-P_i against that complete set
therefore proves absence over the supplied finite base. It does not use a
summation polynomial, a native F5 verdict or a guessed scalar.
"""
from identity import sha256
from oracle import require

POLICY = 'exact-group-three-sum-n17-v1'


class ExactThreeSum:
    def __init__(self, curve, base, summands):
        require(type(summands) is int and summands == 3
                and type(curve.n) is int and 5 <= curve.n <= 17
                and 1 <= len(base) <= 128,
                'exact negative proof adapter domain exceeded')
        # The PDP uses geometric curve points. Final relation coefficients use
        # their cofactor projections; subgroup usability is audited separately.
        # Removing torsion/coset points here would give an incomplete negative
        # proof for the actual solver input.
        require(all(p is not None and curve.decode(list(p)) == p for p in base)
                and len(set(base)) == len(base),
                'exact negative proof base is not distinct finite curve points')
        self.curve, self.base = curve, tuple(base)
        self.pairs = {curve.add(p, q) for p in base for q in base}
        self.rows = []
        self.base_digest = sha256([list(p) for p in base])
        ordered = sorted(self.pairs, key=lambda p: (p is not None, p))
        self.pair_digest = sha256([None if p is None else list(p) for p in ordered])

    def verify_absent(self, target, *, stage, trial, a, b):
        require(target is not None and stage in {'collection', 'target'},
                'invalid exact negative query')
        for point in self.base:
            require(self.curve.add(target, self.curve.neg(point)) not in self.pairs,
                    'false negative: three-point group decomposition exists')
        self.rows.append(dict(stage=stage, trial=trial, a=a, b=b,
                              point=list(target), independently_proved_absent=True))

    def receipt(self):
        return dict(schema_version=1, policy=POLICY, summands=3,
            field_degree=self.curve.n, base_size=len(self.base),
            base_role='exact geometric PDP points before cofactor projection',
            factor_base_sha256=self.base_digest, pair_set_size=len(self.pairs),
            pair_set_sha256=self.pair_digest, repetition_policy='all-indices-with-replacement',
            identity_pair_sums_included=None in self.pairs,
            independently_proved_negative_queries=len(self.rows), rows=list(self.rows),
            scope='postexecution exact finite-group absence; no general F5 soundness or performance claim')
