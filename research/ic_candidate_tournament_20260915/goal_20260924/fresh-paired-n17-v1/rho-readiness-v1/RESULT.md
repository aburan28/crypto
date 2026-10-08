# Disclosed strong-rho readiness control

This is a correctness control on an already exposed public-point workload,
**not** a pre-card source publication or a fresh paired rho result. The
[frozen spec and raw session](session/spec.json) contain one target, one rho
arm, one round, no warmup, and no cross-target collision table. The method is
`rho.signed_frobenius_strong` with 32 lockstep lanes, four distinguished-point
bits and a step-cap factor of 2,000. The target law is
`hash_to_subgroup_v1`, seed 11, index zero. The public point is
`[0x131ff,0x1c6fb]` = `[78335,116475]`; no scalar was supplied to the
solver. This target is excluded from the later fresh card.

The [raw run](session/records.jsonl) recovered scalar **38129** and verified
`[38129]G = Q`. Its [independent ecbench audit](audit.json) checked one
record, reproduced the deterministic run exactly, and found zero problems.
The release binary SHA-256 was
`ead0c3ed756228a339dd4d903bbe90894f1be9be60c5553c0264a4893e9ab96f`.
The session was `ECBS1hf80a0640ca0c`, with workload
`Wfcd1b6a7cdd4` and record `ECR1h882563a8726d`. The recorded online
interval was 108,500 ns, from the first target-dependent walk start through
scalar replay; it is an **L0 ordinary-host diagnostic**. No host-wide CPU
reservation or isolated-host receipt exists, and no IC arm ran on this point
under this control. The rho/IC ratio and speedup are therefore null.

The ecbench curve registry calls this exact representation
`EC1N17Ce1hdfbf24105ef5` and encodes polynomial-bit elements as hex. The
new IC candidate manifests call the same polynomial-bit model
`EC1N17Ckb1hbbe2b5b6b1e6` and encode those integers as decimal. The
explicit crosswalk is field polynomial `x^17+x^3+1` (`0x20009` = 131081),
binary Weierstrass coefficients `[1,1,0,0,1]`, group order 131174,
prime subgroup order 65587, cofactor 2, and generator
`[0xaaad,0x5b2b]` = `[43693,23339]`. Future comparison admission must check
these exact fields and the **same Q**, not infer identity from either short
curve label alone. Input string decoding belongs to setup outside the
one-target online interval for every arm; retain its cost separately.

The run's 1,097.06 group-addition-equivalent count is explicitly a **lower
bound**: canonicalizations, partition hashes, and distinguished-point table
queries/inserts were unpriced. Its `S` must not become a full-operation
comparison without a frozen, matched operation boundary that prices those
terms or labels the omission. This control does not bind a complete source
archive, one-use claim, or post-card execution; those gates remain open.
