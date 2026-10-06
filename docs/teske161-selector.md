# Teske's 161-bit coefficient selection

Source: Edlyn Teske, [*An Elliptic Curve Trapdoor System*](https://doi.org/10.1007/s00145-004-0328-3), *Journal of Cryptology* 19 (2006), §§2.1, 4 and Appendix A.

This native selector reproduces **the first structural step** of the published construction over
(K=mathbb F_{2^{161}}=mathbb F_2[z]/(z^{161}+z^{18}+1)).
It does not certify a usable trapdoor or recover a logarithm.

## Run

```sh
cargo build --release --bin teske161_select
./target/release/teske161_select --component w1 --genus 7 --seed 20261006
./target/release/teske161_select --component w2 --genus 8 --seed 20261006 --a 1
cargo test --bin teske161_select
```

Output is JSON with the exact polynomial-basis coefficient (b), the seed, component, augmented magic number, trace and genus branch. The seed is for reproducibility, not private key generation.

## Selection certificate

Let (sigma(c)=c^{2^{23}}). Its order is seven. The three factors of
(x^7+1=(x+1)(x^3+x^2+1)(x^3+x+1)) define (W_0,W_1,W_2).
For uniformly sampled (uin K), the selector applies one of these linear maps:

| Chosen component | Projection polynomial in (sigma) | Kernel equation checked |
| --- | --- | --- |
| (W_1) | (1+sigma^2+sigma^3+sigma^4) | (1+sigma^2+sigma^3) |
| (W_2) | (1+sigma+sigma^2+sigma^4) | (1+sigma+sigma^3) |

These maps are surjective onto the respective 69-dimensional spaces. The selector rejects zero. For genus 8, it also samples a nonzero (W_0) component by the relative trace (sum_{i=0}^{6}sigma^i(u)); genus 7 uses (W_0=0). It then checks that (b) is in (W_0oplus(W_i\setminus\{0\})), computes Teske's actual magic number as the (mathbb F_2)-rank of the **pairs** ((1,\sqrt{\sigma^i(b)})), and checks the relative trace of (b). Theorem 1 makes trace zero equivalent to genus 7 for these candidates; otherwise the genus is 8.

**Existing generic-helper caveat:** `ec_trapdoor::magic_number_full` ranks only the square-root orbit (plus a separate coefficient contribution). Its `ghs_genus_with_type` uses the orbit length as its type test. Those are not Teske's augmented-pair rank and relative-trace branch, so `ghs_screen` can report a different genus for the very same 161-bit coefficient. Use the selector's paper-specific checks for this construction; reconcile the general helper with a separate mathematical audit before treating its genus as exact here.

## What remains for the full construction

1. Count both twists, prove a large prime-order subgroup with cofactor 2 or 4, and check the embedding degree.
2. Compute the Frobenius discriminant and certify Teske's squarefree, magnitude and class-group conditions.
3. Walk split odd-prime isogenies from the chosen source and retain every certified kernel and signed point map. The existing `binary_isogeny` code walks (j)-invariants; `binary_velu` provides explicit maps in its bounded toy-field implementation. Neither is an end-to-end 161-bit route certificate.
4. Construct the smooth genus-7/8 curve and an explicit subgroup-preserving GHS map, then account for descent, index calculus and all setup against a matched rho reference. The existing `ghs_descent` code executes only its (m=1) case end to end.

A structural candidate alone says nothing about a NIST curve. This field has a degree-23 intermediate field and this particular construction does not transfer to arbitrary prime-degree binary fields or prime-field curves.
