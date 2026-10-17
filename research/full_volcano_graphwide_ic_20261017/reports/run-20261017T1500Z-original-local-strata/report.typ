#set page(paper:"a4",margin:13mm)
#set text(size:9pt)
= Original A--E: local End(E) conductor certificates
17 October 2026. Input certification only; no new ECDLP timing run.
#image("coverage.svg",width:100%)
#table(columns:(.3fr,1.8fr,1.2fr,.6fr,.6fr,.6fr),inset:3pt,
 [Case], [Prime field p], [Source conductor], [Prefix vertices], [Exact vertices], [Observed strata],
[A], [126100789566376697], [268435456], [1538], [512], [10],
[B], [45862249674083071], [129140163], [1616], [729], [7],
[C], [1637958981148213249], [967458816], [1869], [517], [24],
[D], [814409292524040673], [544195584], [1881], [509], [24],
[E], [1975950516940704721], [644972544], [1881], [509], [24],
)

The conductor is the requested index of End(E) in the maximal CM order.
These source values are locally determined from complete rational-kernel
counts and the standard ordinary-volcano floor-distance correspondence;
they are not merely copied Frobenius conductors. For each prime, distance
to known floors bounds height above; distance to known floors or unexpanded
frontier bounds it below. Only equal finite bounds produce exact labels.
Unknown vertices are never declared floors or omitted from coverage.

The same calculation passed all available complete-truth bound/label checks
on the three complete additional graphs at visibility cutoffs 0, 1, 17, 1024
and all. The generic certify mode reproduces the earlier C certificate byte
for byte. It shares isogeny arithmetic and the standard theorem; this is not
an independent arithmetic implementation or a new mathematical theorem.

All five frozen prefix DBs are read immutably. Every output binds its DB
SHA-256 and existing map-replay receipt. The original databases and global
stratum status are unchanged. Actual local certificate rows and per-conductor
counts are in the accompanying archive. No prefix is labeled a complete class.
Original full A--E total 7,833,658,967 expected vertices and 515 strata;
historical B has a separate complete delivery. A/C/D/E remain incomplete.

The largest prime subgroup and embedding degree are separate properties.
In particular C has largest prime order 7 despite its 61-bit field and
source conductor 967,458,816. These local labels do not establish a speed
advantage, a vulnerability caused by conductor, a full benchmark matrix or
finite-field solver support for all models. No deployed-security claim.
