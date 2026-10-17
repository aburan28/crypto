# Local dashboard browser limitation

The static dashboard generator and evidence-pin check pass. The local Chromium
render check cannot navigate to the repository's `file://` page: the managed
browser returns `chrome-error://chromewebdata/` with “This page is blocked —
Your organization doesn't allow you to view this site.” Its zero-chart assertion
therefore does not describe the application DOM.

The first launch also encountered this container's unavailable SUID sandbox;
the normal container fallback reached the managed browser policy block. Both
diagnostics are retained. The policy was not disabled or circumvented. The
existing GitHub dashboard workflow performs the actual desktop/mobile render
checks in its normal runner environment after the generated-page repair.

This is a local validation limitation, not a passing browser-render receipt.
