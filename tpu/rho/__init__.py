"""
tpu/rho — batched Pollard rho (parallel collision search) for prime-field
ECDLP, with the modular multiply routed through int8 matmuls.

Provenance: Python (JAX) by explicit user direction for the tpu/ backend
(2026-10-01), overriding AGENTS.md's no-Python rule for this work.
Correctness-only; no device has run it and no speed claim is made.
"""
