# Result publication integration, 4 October 2026

The original F5 registration remains consumed and closed. Integration with
main `a9c8058b945d0f4b20daf6bd41ebab6fa9ac581b` changes no file under
`result-v1`, no original source seal, execution receipt or frozen audit, and
starts no scientific worker or auditor. The result remains a disclosed n17
correctness control, with no new natural yield, fresh qualification or speedup.

The previous PR head `32dab3e36f81a77664e2215a607f3bc72637f8e8` passed all 49
applicable CI checks, including native result replay on Linux and macOS. Main
advanced while those checks ran. The three conflicts were the presentation
renderer, generated dashboard and generated panel index. The resolution keeps
the source-bound native F5 result pin and diagnostics, main's Lab browser link,
all newer source-pinned panels and every historical panel. Regeneration checks
pass, with 261 navigation panels. The original historical ledger and regime
summary remain preserved.

Local integrated validation passed all 47 website structure/navigation tests.
The real Chrome check passed at 1280, 390 and 320 pixels, in dark mode and with
JavaScript disabled, with zero script errors and external requests. The checked
page SHA-256 is
`ac99d3a7439485d4abd7f7c755b63623cca26fc205a9801cd68e0499e0504747`.

Two operational failures preceded those checks: test discovery was first
pointed at the nonexistent `tests/site` directory; the corrected discovery
uses `scripts/site/test_build.py`. The sandboxed Chrome launch then exceeded
its existing launch deadline. A fresh-profile run with normal local browser
process permissions passed. No timeout, assertion or scientific gate was
relaxed. Compact raw logs are retained beside this note. The combined website
and sandbox-browser log is gzip-compressed to preserve its exact original
whitespace; `gzip -cd integration-site-and-sandbox-browser-20261004.txt.gz`
extracts it without modifying the evidence.

The integrated PR head still needs its own applicable CI and review gates.
Success of the earlier head is not acceptance of this integration. Browser and
correctness checks are not CPU timing calibration. Under the updated isolation
rule, no controlled CPU wall-time speedup is admitted without the required
auditable host-level receipt; this control's speedup remains null.
