"""ICMS: the index-calculus measurement standard, version 1.

The standard is docs/ic/measurement/README.md.  This package implements it:

* spec.py        load and validate an experiment specification (YAML or JSON)
* canonical.py   canonical JSON, hashing, and every identifier the standard mints
* environment.py the host capsule: stable facts (hashed into env_class_id)
                 and volatile facts (sampled around every run)
* execute.py     one pinned, sampled execution of a measured command
* gates.py       the isolation gate: which level (L0-L3) a run earned
* registry.py    controlled vocabularies: units, references, windows, phases
* adapters/      compile a spec into a producer's command line and parse its
                 report into the record's metrics
* session.py     a measurement session: lock, preflight, eviction,
                 interleaved repetitions, records
* compare.py     admission rules and paired ratios between two run sets

Standard library plus PyYAML (only to read YAML specs; JSON specs need nothing).
"""

STANDARD = "icms/v1"
SPEC_SCHEMA = "icms.spec/v1"
RECORD_SCHEMA = "icms.record/v1"
CAPSULE_SCHEMA = "icms.capsule/v1"
SESSION_SCHEMA = "icms.session/v1"
COMPARISON_SCHEMA = "icms.comparison/v1"
