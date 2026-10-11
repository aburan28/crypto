"""Experiment 3 (kill-or-keep): an attacker who tries to recover the secret
endomorphism from the public eigenbasis guesses Delta' (hence s0 = sqrt(-Delta')
mod 2^(e+2) and the Kani scalars c, c' with c^2 + c'^2 = 2^e - d^2 Delta').
Every such guess yields a maximal isotropic kernel, so the chain runs.  Does a
wrong guess fail only at the end (no early-abort signal) or does the failure
step correlate with the number of correct low-order bits of Delta'?"""
import sys, random, time
from sage.all import ZZ, set_random_seed, randint, two_squares
import eigen_seta_toy as T
from kani_chain import chain_evaluate

def consistent_guess(Dg, N2, ep, d):
    try:
        cc, cp = two_squares(N2 - d * d * Dg)
    except ValueError:
        return None
    if cc % 2 == 1:
        cc, cp = cp, cc
    if cp % 2 == 0:
        return None
    s0 = T.sqrt_mod_2power(-Dg, ep + 2) % N2
    return ZZ(cc), ZZ(cp), s0

def run(e, f, k, variant, n_rand=40, seed=1):
    set_random_seed(seed); random.seed(seed)
    p, c = T.find_prime(e, f, k)
    pk, sk, par = T.keygen(p, e, f, k, variant=variant, verbose=False)
    N2, d, ep = par["N2"], par["d"], par["e"]
    ct, wit = T.encaps(pk, par)
    EC, S, Tt = ct["EC"], ct["S"], ct["T"]
    Dp = sk["Dp"]
    print(f"{variant}: p {p.nbits()} bits, chain length {ep}, true Delta' = {Dp}")
    X1, _ = T.torsion_basis(EC, d, (p + 1) // d)
    O = EC(0)
    guesses = []
    # guesses agreeing with Delta' on the low j bits, for a range of j
    for j in range(3, ep - 1, max(1, ep // 7)):
        tries = 0
        while tries < 200:
            tries += 1
            Dg = (Dp % 2 ** j) + 2 ** j * randint(1, max(1, (N2 // (d * d)) // 2 ** j))
            if Dg == Dp or Dg % 8 != 7 or d * d * Dg >= N2:
                continue
            g = consistent_guess(Dg, N2, ep, d)
            if g:
                guesses.append((j, Dg, g)); break
    for _ in range(n_rand):
        while True:
            Dg = 8 * randint(0, (N2 // (d * d)) // 8) + 7
            if Dg != Dp and d * d * Dg < N2:
                g = consistent_guess(Dg, N2, ep, d)
                if g:
                    guesses.append((0, Dg, g)); break
    guesses.append(("true", Dp, (sk["c"], sk["cp"], sk["s0"])))
    for j, Dg, (cc, cp, s0) in guesses:
        lamS = (d * s0 + cc) % N2; lamT = (cc - d * s0) % N2
        K1 = (cp * S, (-lamS) * S); K2 = (cp * Tt, (-lamT) * Tt)
        assert S.weil_pairing(Tt, N2) ** (cp * cp + lamS * lamT) == 1
        t0 = time.time()
        try:
            (E1, E2), imgs, log = chain_evaluate(EC, EC, K1, K2, ep, [(X1, O)])
            ok = E1.is_isomorphic(EC) and E2.is_isomorphic(EC)
            verdict = f"ran {len(log)}/{ep} steps, ended on product, codomain is E_C x E_C: {ok}"
        except RuntimeError as ex:
            verdict = "ran all steps, ended on a Jacobian (no final split)"
        except Exception as ex:
            verdict = f"failed: {type(ex).__name__}"
        tag = "TRUE" if j == "true" else (f"low {j} bits right" if j else "random")
        print(f"  Delta' guess {tag:>18}: {verdict} ({time.time()-t0:.2f}s)")

if __name__ == "__main__":
    e = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    f = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    variant = sys.argv[4] if len(sys.argv) > 4 else "cm"
    run(e, f, k, variant)
