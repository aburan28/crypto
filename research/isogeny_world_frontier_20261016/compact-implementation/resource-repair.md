Native resource-admission repair
===============================

The public supervisor required positive RSS before acknowledging a waiting
compact worker, but previously accepted a later Some(0) sample as usable
telemetry. The changed predicate rejects zero during the entire supervised
work phase. Missing/zero readings use the existing 100 ms fail-closed grace;
the worker's done marker still bounds the completed mathematical phase.
There is no demonstrated false mathematical acceptance from this gap.

The internal sampler dependency permits deterministic native controls without
adding any CLI or environment override to the production monitor. On Unix,
a real shell fixture waits for the same acknowledgement, records entry, then
remains alive until the supervisor stops it. It performs no cryptographic work.
Four cases inject missing and zero RSS before and after acknowledgement. Each
returns RESOURCE_MONITOR_UNAVAILABLE; pre-ack cases create no acknowledgement
or entry, post-ack cases retain both. fixture.json explicitly marks the RSS
values synthetic: they are not host-memory measurements.

The focused run passed (exit 0), with operational durations 102/121/161/158 ms.
Its raw records are under resource-monitor-retake-v1/. A separate full native
suite passed 91 tests, zero failures: 3 field checks, 82 main-binary tests,
4 registry-model tests and two CLI integrations. Its distinct four telemetry
controls are under resource-monitor-full-v1/, and positive/forged-certificate
CLI records are under native-controls-resource-retake-v1/. Both logs and exits
are retained. These are correctness controls; no comparative speedup is claimed.
The injected-process controls are Unix-only; Windows is not claimed as executed
here. Existing public CLI tests remain portable and await their current CI head.

The precommit binary's source/dirty fields and hashes are preserved beside
this note. Historical clean build9a2 and all its raw controls remain unchanged;
it is superseded for the proposed M8191 scientific launch by this resource
repair. A new committed clean build, source-bound review/control retake, custody
and frozen protocol are required. No M8191 target has been executed.

The report and controls diagram were refreshed; expanded-map challenge scores,
season, canonical IC/rho graphs and the curve registry are unchanged. Root src
and Cargo files are unchanged; the already-passed root release-library gate
(4148 passed,117 ignored,zero failures) remains the relevant root validation.
