"""Fiber invisibility check (RESEARCH_ISOGENY_FIBERS_ECDLP.md, section 13):
every F_p-rational isogeny of degree coprime to r has image containing the
prime-order subgroup G, the intersection of the images is G, and the
connecting class delta(Q) in E2(F_p)/phi(E1(F_p)) is constant on cosets of G.
Toy curve of section 5: y^2 = x^3 + x + 3 over F_211, N = 228 = 12 * 19.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fiber_probe import Ext, Curve, velu, tonelli

p = 211; F = Ext(p, 1, [0, 1]); E1 = Curve(F, F.emb(1), F.emb(3))
pts = E1.points_prime(); N = len(pts); r = 19; h = N // r
G = {str(E1.mul(h, P)) for P in pts}                     # [h]E(F_p) = G
assert len(G) == r
T2 = next(P for P in pts if P is not None and F.is_zero(P[1]))
T3 = next(P for P in pts if P is not None and P[0][0] == 77)
out = {}
for name, ker in (("2-isogeny", [T2]), ("3-isogeny", [T3])):
    E2, phi = velu(E1, ker)
    image = {str(phi(P)) for P in pts}
    G2 = {str(E2.mul(h, Q)) for Q in E2.points_prime()}
    pts2 = E2.points_prime()
    # delta(Q): coset of the image in E2(F_p); label by a canonical representative
    cosets = {}
    for Q in pts2:
        key = str(Q)
        if key in cosets: continue
        label = len(cosets)
        # coset Q + image
        for P in pts:
            cosets[str(E2.add(Q, phi(P)))] = label
    # constancy on cosets of G2 = [h]E2(F_p) (the prime-order subgroup of E2)
    Gpts = [Q for Q in pts2 if str(Q) in G2]
    constant = all(len({cosets[str(E2.add(Q, g))] for g in Gpts}) == 1 for Q in pts2)
    out[name] = dict(index=len(pts2)//len(image), G_in_image=G2 <= image,
                     delta_constant_on_G_cosets=constant, n_cosets=len(set(cosets.values())))
# Intersection of images on E1: for each rational isogeny phi: E1 -> E2 the dual
# phi_hat: E2 -> E1 has image [deg phi] E1(F_p) (phi_hat o phi = [deg]); the
# intersection over all rational isogenies of degree coprime to r is [h]E1(F_p) = G.
images_on_E1 = []
for ker in ([T2], [T3]):
    deg = 2 if F.is_zero(ker[0][1]) else 3
    img = {str(E1.mul(deg, P)) for P in pts}
    images_on_E1.append(img)
inter = set.intersection(*images_on_E1)
out["intersection_of_dual_images_on_E1"] = dict(size=len(inter), contains_G=G <= inter, equals_G_after_h=(G == {str(E1.mul(h, P)) for P in pts}))
print(json.dumps(out, indent=1))
