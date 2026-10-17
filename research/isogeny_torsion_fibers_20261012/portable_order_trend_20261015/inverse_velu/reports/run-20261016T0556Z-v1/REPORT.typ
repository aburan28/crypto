#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Inverse-Vélu V1 splitter: failed smoke child
*Public cell `p=1237,r=1201`; round zero, one target.*

The V1 schedule used 32 base-field shifts `[i,0]` to split the
rational factor of the inverse cubic by quadratic character. The
native child hit its frozen cap and exited 101 before relation
collection or scalar recovery. Its first stderr, unchanged-binary
retry, exit file, and source/binary SHA-256 are retained beside the
versioned protocol. The first shell wrapper did not record its exit
because `status` is read-only in zsh; the retry did.

#image("DIAGRAM.svg", width: 100%)

For roots and shifts in the base field, every nonzero shifted value
is a square in the quadratic extension. That makes this character
test unable to split such a cubic. V2 changes only the shift family
to `[i,1]`; the V1 failure is not treated as a missing fiber or an
ECDLP cost result. See `REPORT.md` and `../../PROTOCOL_V2.md`.
