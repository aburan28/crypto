# Stage 199 generic-receipt correction

The first charged finalizer attempt failed after the unset solver replay had
already passed. The generic development-receipt audit incorrectly required
every receipt to have the solver's zero exit, no-timeout status, and 360-second
watchdog. That rejected the legitimate 30-second layout receipt and would also
exclude failed attempts from additive accounting.

The correction keeps schema, artifact, metrics, resource, and positive-watchdog
checks in generic receipt replay. Solver-specific zero-exit, no-timeout, and
360-second requirements now live in the unset solver-command validator. The
failed finalizer receipt is preserved and charged. No solver measurement,
selected policy, input, or decision changes.

