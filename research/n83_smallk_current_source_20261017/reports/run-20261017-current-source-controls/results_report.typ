#set page(width: 297mm, height: 210mm, margin: (x: 13mm, y: 11mm))
#set text(font: "Helvetica Neue", size: 9pt)
#set heading(numbering: none)

= N83 current-source factor-base comparison

K1280/root11 dynamic reached 4.73325 ms online and 146.960 s receipt wall
on the previously exposed public point. The same-source K2543/root11 control
used 18.983375 ms online and 401.115 s receipt wall. K1280's observed online
interval was 4.01 times shorter on this point.

#figure(image("controls.svg", width: 100%), caption: [Verified one-target online intervals and cold receipt walls. One process per arm on a shared macOS host. Values are diagnostics on one exposed point; uncertainty is not estimated.])

The eight rows use the same IC source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and local release binary SHA-256
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
All full-base, four-point relation, scalar, phase and candidate-identity
replays passed. The table and plot are derived from the versioned receipt
and replay files listed in `plot_data.json`.

The next selection step is a fresh-point Tier 1 alternating same-point
IC–rho comparison. The canonical IC/rho dashboard remains unchanged by
this local control.
