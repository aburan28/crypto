#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)

= Selected order-1,201 rational-fiber audit

The frozen ascending-prime selector chose $(p,r,h)=(1237,1201,1275)$.
The Rust runner enumerated every source point and destination image;
an independent Rust implementation reconstructed every fiber and all
$1201^2 = 1,442,401$ ordered subgroup pair identities. Both processes
exited 0. Their source/image counts are 1,531,275 and 510,425,
respectively, with three preimages at every rational image point.
The identity fiber is $\{O,T,-T\}$ and the order-1,201 subgroup maps
injectively. The runner/checker peak RSS values were 199,163,904 and
367,984,640 bytes, below the 512 MiB cap.

#image("FIBER_V1.svg", width: 100%)

Exact commands, source/binary and receipt hashes, stderr, the full-map
and subgroup digests, and links to the native sources are in
`FIBER_RUN.md`, `FIBER_1237_1201_V1.json` and its `.check.json`.
This proves the selected finite fiber and map-law identities only.
The cold same-target solver, complete portable cost and isolated
performance gates follow separately.
