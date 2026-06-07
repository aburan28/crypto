"""Shared cheap-invariant scoring (Phase 4 features) for candidate curves."""

from ecfp import Curve


def cheap_score(E: Curve):
    p = E.p
    j = E.j_invariant()
    if j == 0:
        aut = 6
    elif j == 1728 % p:
        aut = 4
    else:
        aut = 2
    return {
        "j_is_0": j == 0,
        "j_is_1728": j == 1728 % p,
        "a_is_minus3": E.a == (-3) % p,
        "a_hamming_weight": bin(E.a).count("1"),
        "b_hamming_weight": bin(E.b).count("1"),
        "automorphism_group_size": aut,
        "special": (j in (0, 1728 % p)) or E.a == (-3) % p,
    }
