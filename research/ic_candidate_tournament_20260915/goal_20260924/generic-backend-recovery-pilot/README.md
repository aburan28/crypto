# F4/F5/SAT disclosed recovery pilot

Status: **all twenty registered jobs completed or reached their fixed process
cap**. The [results](RESULTS.md) and [complete evidence bundle](evidence.tar.gz)
retain one verified n17a1 F5 recovery, four process timeouts and fifteen
bounded-incomplete reports. The combined F4/F5-and-SAT admission gate failed:
neither n17a1 SAT arm returned a report. The [protocol](PROTOCOL.md) and
[panel](panel.json) froze one source, five previously disclosed public points,
four solver backends and per-cell budgets before execution. The
[runner](../../run_generic_backend_recovery_pilot.py) executed every job once;
the [evidence replay test](../../test_generic_backend_recovery_evidence.py)
rechecks the retained process and scientific admission receipts. No job was
retried.

This is the next feasibility gate after the
[negative one-query pilot](../generic-backend-disclosed-pilot/RESULTS.md).
The verified complete solve is on a disclosed point and does not establish
fresh-target yield, a matched rho ratio or tournament promotion.
The separate v2 seed `2026092902` campaign remains the sole permitted
measurement of its registered workload and must not be dispatched again.
