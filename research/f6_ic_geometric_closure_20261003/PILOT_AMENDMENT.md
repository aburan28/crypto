# Pilot amendment before timing

The registered protocol proposed three arms: inherited F4, F6-IC and rho.
After the protocol was committed, the repository's current `AGENTS.md` §12
was checked against this older tournament worker. It requires new comparisons
of ECDLP method classes to use `ecbench`; the prepared worker does not yet
have an `ecbench` adapter. No F6-IC pilot timing has been run as of this
amendment.

The bounded pilot therefore compares **only the two IC PDP variants** through
the existing prepared worker. It keeps the frozen n17 mathematical state,
target seed `20261003019`, algorithm seed `20261003034`, one worker, the same
public point, limits, five exclusive phases, source/binary/input freeze, and
600-second cap. It runs six fresh processes in this fixed order:

1. inherited F4, F6-IC;
2. F6-IC, inherited F4;
3. inherited F4, F6-IC.

The full one-target online intervals and correctness records are useful
candidate diagnostics. The local Mac has no host-level CPU isolation receipt,
so even a favorable paired time remains **exploratory**. No F6-versus-rho
measurement or speedup is claimed by this pilot. A cross-method comparison
must add the candidate to `ecbench`, use the same public target and resource
envelope as strong rho, produce a sealed replayed session, and satisfy the
independent and isolated-host gates before any controlled wall-time claim.

All original stop rules still apply. The six-process completion gate replaces
the original nine-process gate, and no target or limits change after results.
