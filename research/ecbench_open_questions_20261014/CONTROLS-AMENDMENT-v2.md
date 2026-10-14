# Structural-control amendment 2

2026-10-14, before the second structural-control run. Measurement specs and
the required E0 degree83 gate are unchanged. The first structural attempt
failed in the wide E1 constructor with `#E / 2 = 4835703278455440964278031
is not a usable prime`; its log is retained as `controls-v1-failed.log`.
No completed control receipt was written, so that attempt establishes no
complete structural panel.

Use the exact existing E1 degree83 catalog representation
`EC1N83Ce1h35c7ced6b497`, model `icv1-f2m83-t6151469093347-9eba6cd5`,
modulus z83+z7+z4+z2+1, r=8569786107849059, cofactor1128547018,
G=(0x5789eebc7f493ecfe93af,0x755ed5c31e037ed90b549). The original
question asked for the registered E1 model; it did not authorize asserting
cofactor2. Preserve the constructor refusal and identify this exact subgroup
in the successful pullback receipt. It is a smaller subgroup than E0's N83
gate and never substitutes for that gate.

Preserve each structural stage outcome independently. A failure in one stage
must not discard unrelated observations. The second receipt has schema
`ecbench-ten-questions-controls/v2` and names this amendment. All original
acceptance gates, including independent source-order proof and absence of a
general Jacobian norm/DLP solver, remain visible.
