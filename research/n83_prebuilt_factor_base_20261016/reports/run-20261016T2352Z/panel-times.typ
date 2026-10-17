#set page(width: 170mm, height: 105mm, margin: 6mm)
#set text(size: 8pt)
#let bar(seconds, paint) = rect(width: seconds * 0.52mm, height: 5mm, fill: paint, radius: 1mm)

*N83 retained factor-base panel: three local repetitions*  \
Seconds, local L0 host; each panel contains 54 compressed objects (83,328,324 bytes).

#table(
  columns: (25mm, 98mm, auto),
  inset: 4pt,
  align: (left, left, right),
  [*Interval*], [*0 → 180 s; bar length = 0.52 mm/s*], [*Seconds*],
  [Build 1], [#bar(59.517962583, rgb("2563eb"))], [59.518],
  [Build 2], [#bar(74.022112708, rgb("2563eb"))], [74.022],
  [Build 3], [#bar(169.554772250, rgb("2563eb"))], [169.555],
  [Replay 1], [#bar(46.66, rgb("16a34a"))], [46.66],
  [Replay 2], [#bar(52.15, rgb("16a34a"))], [52.15],
  [Replay 3], [#bar(31.43, rgb("16a34a"))], [31.43],
)

Blue: native construction + compression + object writes (internal interval).  \
Green: separate full generic arithmetic replay of a saved manifest against
byte-identical reconstructed objects (outer interval).  \
All six producer/replay data rows completed; the original `/usr/bin/time -l`
wrapper failed after each build because sandboxed `sysctl kern.clockrate` was
unavailable. Build and replay perform different work; no speed ratio is inferred.

*Run ID:* `run-20261016T2352Z`  \
*Source:* `runs/run-{1,2,3}/manifest.json` and `replay-recovery.stderr`.
