"""TPU (JAX/Pallas) index-calculus backend.

Two stages, each expressed in the bit-matrix-multiply-mod-2 shape a TPU's
systolic array executes natively:

* :mod:`tpu.ic.pairtable` -- stage 1, pair-table relation collection.
* :mod:`tpu.ic.linalg`    -- stage 2, GF(2) relation algebra.

:mod:`tpu.ic.reference` is the scalar oracle every kernel is checked
against.  See ``tpu/README.md`` and ``tpu/protocol/RESEARCH_TPU_IC.md`` for
the honest fit analysis and the (all-pending) stage-diagnostic table.

This directory is Python/JAX by explicit user direction, which overrides
``AGENTS.md``'s no-Python rule for this work; it is labelled as such so its
provenance is never in doubt.
"""

from . import curve, field, instances, linalg, pairtable, reference  # noqa: F401

__all__ = ["curve", "field", "instances", "linalg", "pairtable", "reference"]
