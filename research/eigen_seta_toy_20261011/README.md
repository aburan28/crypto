# Eigen-Séta toy implementation status

The source is `isogeny_algos/src/eigen_seta_toy.rs`. Run
`cargo run --release --example eigen_seta_toy` from `isogeny_algos/`.

This is a bounded correctness experiment on the public CM curve of
discriminant -7. It is **not** an implementation of the full proposal and
must not protect data. Its public sender-index search takes at most 81
isogenies, so anyone can derive the toy shared key.

| Requested part | Status | Evidence and gap |
| --- | --- | --- |
| Key generation | Toy complete | Computes a genuine norm-2 CM endomorphism and its 2-power eigenbasis; the CM curve and its endomorphism ring are public. Random-curve key generation is absent. |
| Encapsulation | Toy complete | Builds the 3-power Vélu isogeny, applies deterministic odd diagonal masks, and hashes a canonical ciphertext transcript. |
| Decapsulation | Toy complete by public exhaustive recovery | Re-encrypts every sender index. The secret endomorphism gives no decryption advantage in this implementation. |
| First Kani half | Toy complete on the tested fixture | Constructs the ciphertext graph, uses an explicit 4-division table, handles product-to-product steps, then runs the existing theta chain. |
| Second Kani half and component maps | Not implemented | The existing theta chain returns an intermediate Jacobian. It does not construct the dual chain, identify the middle Jacobians, or recover the kernel of the dual sender isogeny. |
| Full-size parameters and security | Not established | The toy uses a 16,127 prime and a tiny sender degree. No constant-time arithmetic, CCA transform, random-curve keygen, or security reduction exists here. |

The design note derives independent sender index `s` and masks `a,b` from a
seed, then says to recover the seed after recovering `s`. That inversion is
not specified. This toy instead uses `s` itself as the seed and derives the
masks from it. Re-encryption is consequently well-defined, but public
enumeration of `s` is trivial at these parameters.

For the tested `p=16127`, `N=2^6`, `d=3^2` fixture, the first Kani half
needs an initial diagonal quotient followed by several factorwise 2-isogeny
quotients before the existing theta chain applies. Its final quotient is a
Jacobian, matching its role as the *first* half. The second half remains an
open implementation requirement; passing the KEM round trip does not fill it.
