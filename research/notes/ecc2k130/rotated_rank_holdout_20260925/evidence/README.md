# Frozen rank and point-only holdout evidence

Preregistered source/input head:
`c1ae12eec1469241da560cd329184e1843476915`.
`FROZEN.json` SHA-256:
`bcf0ff3fe38c50fcfcece354a1eb3ec09e8ac5fa021c98e55446e4872f97de4a`.
Immutable #766 row archive SHA-256:
`863f50ea872ce01b1999d0ef67bc0341fa0c764f2ad4187d16736768b0ab9b3d`.
The frozen point-only Q file and sealed k file are separate; neither was
changed after the preregistration commit. The selection transcript preserves
all SHA-ordered candidates scanned before 16 positive Q per arm.

Accepted `raw.tar.gz`: 4,134 bytes, SHA-256
`3b483e486eb1dec9d34a81653cb6b568b74ada33e8c05305304810ec87ba3c6a`.
`receipt.json`: SHA-256
`b7fd19ad37caaf58bb148b1c1bde5f822483be1da6b00d7f6e27feb82836c08b`.
The archive contains modular training counts and base logs, the complete
charged archive-oracle response, point-only recovered scalars and the
independent sealed-label replay report. The receipt records every child
argv, wall/CPU/RSS diagnostics, all raw-file hashes, and partial-file hashes
on failure. The source/input freeze includes the input-construction receipt.

Run `python3 ci_replay.py` from this directory. It validates the protocol-
pinned source/input hashes, regenerates the SHA split, checks the output
archive and receipt, then independently recomputes rank, point-log recovery
and sealed-label comparisons using a separate bit-serial/Fermat group law.
The archive-only path never reruns the rank producer or point-only child.
