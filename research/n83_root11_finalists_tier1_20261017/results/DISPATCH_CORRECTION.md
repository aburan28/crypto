# Provider-censoring correction before cloud execution

The first sealed dispatch at commit `f56cbf870eed3a75861551826df6a9304a781bbe`
had dispatch SHA-256
`77cfee4e81259f45ad4fc00698b41603e7b22dce046cb1e481f901077cedcd59`.
GitHub Actions run `39060689296` was cancelled while all four matrix jobs were
queued. The provider reports all four job conclusions as `cancelled`; none
entered the producer phase.

The provider identity binder originally required every job to conclude
`success`. A calibrated but noisy stratum writes a replayable
`BOUNDED_UNKNOWN` analysis, exits 2 and causes a provider job conclusion of
`failure`. The successor binder admits that job only as a censored stratum
when its archived analysis matches independent local replay, its local exit
is 2 and its analysis status is `BOUNDED_UNKNOWN`. The full-panel selector
still requires all strata admitted and cannot select a candidate from this
path. The frozen targets, exact candidates, schedule and statistical rule
are unchanged. The original commit and cancelled run retain the earlier
dispatch.
