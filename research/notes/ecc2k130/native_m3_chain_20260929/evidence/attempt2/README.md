# Expanded run before canonical wrong-target correction

This complete second held run used source commit
`d068fcffb2725a06b2f148c046983a84511859d2` and its copied
`FROZEN.json`; all toy and leaf checks passed, and the source-specific
independent replay passed. It checked all 27 triples and found a fully lifted
generic/generic SAT model.

Final diff review found that the raw wrong-target selector iterated the
oracle's `[O, finite points...]` list, whereas the preregistered protocol
specified canonical lexicographic `(o,x,y)` order, where finite points precede
`O`. This affects the identity of the 27 wrong targets, even though every
archived wrong target was rejected. The raw CNFs, solver streams, rows,
receipts and hashes are preserved here. A new source lock and a fresh final
run implement the protocol's exact ordering. This archive is not used for the
final semantic decision or claimed as natural-target yield.
