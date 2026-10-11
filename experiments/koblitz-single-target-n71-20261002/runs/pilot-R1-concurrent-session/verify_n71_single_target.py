from pathlib import Path
import json
from sage.all import EllipticCurve, GF, Integer, PolynomialRing, factor

run_dir = Path("/Volumes/SSD990/crypto/experiments/koblitz-single-target-n71-20261002/runs/R1")
ic = json.loads((run_dir / "ic.jsonl").read_text().splitlines()[0])
rho = json.loads((run_dir / "rho.jsonl").read_text().splitlines()[0])
fixture = json.loads((run_dir / "public_target.json").read_text())
fixture_scalar = Integer((run_dir / "target_scalar_validation_only.txt").read_text().strip())

R = PolynomialRing(GF(2), "z")
z = R.gen()
modulus = z**71 + z**5 + z**3 + z + 1
assert modulus.is_irreducible()
F = GF(2**71, name="u", modulus=modulus)
u = F.gen()

def decode(value):
    value = Integer(value)
    out = F(0)
    for i in range(71):
        if (value >> i) & 1:
            out += u**i
    return out

E = EllipticCurve(F, [1, 0, 0, 0, 1])
G = E(decode(fixture["generator"][0]), decode(fixture["generator"][1]))
Q = E(decode(fixture["public_target"][0]), decode(fixture["public_target"][1]))
O = E(0)
r = Integer(fixture["subgroup_order"])
lam = Integer(5336382444749771)
curve_order = Integer(2361183241386169526132)
cofactor = Integer(428276)
trace = Integer(48653080717)
assert curve_order == r * cofactor == 2**71 + 1 - trace
frobenius_G = E(G[0]**2, G[1]**2)
assert lam * G == frobenius_G
assert pow(int(lam), 71, int(r)) == 1 and lam % r != 1

assert G != O and Q != O
assert r * G == O
for prime, exponent in factor(r):
    assert (r // prime) * G != O, "generator order is smaller than claimed subgroup order"
assert Integer(ic["recovered_scalar"]) == fixture_scalar
assert Integer(rho["recovered_fixture_scalar"]) == fixture_scalar
assert ic["published_q"] == fixture["public_target"]
assert rho["published_q"] == fixture["public_target"]
assert ic["generator"] == fixture["generator"]
assert rho["generator"] == fixture["generator"]
assert Integer(ic["recovered_scalar"]) * G == Q
assert Integer(rho["recovered_fixture_scalar"]) * G == Q
assert rho["verified"] is True and rho["reference_group_validation"] is True
assert ic["group_verified"] is True and ic["exit_code"] == 0

result = {
    "kind": "sage_independent_single_target_scalar_replay",
    "field_degree": 71,
    "field_polynomial": "x^71 + x^5 + x^3 + x + 1",
    "curve": "y^2 + x*y = x^3 + 1",
    "subgroup_order": str(r),
    "curve_order_record_consistency_verified": True,
    "curve_trace": str(trace),
    "frobenius_eigenvalue_mod_r": str(lam),
    "frobenius_action_on_generator_verified": True,
    "frobenius_order_on_subgroup": 71,
    "signed_frobenius_orbit_size": 142,
    "subgroup_order_factorization": str(factor(r)),
    "target_count": 1,
    "public_target": fixture["public_target"],
    "fixture_scalar_validation_only": str(fixture_scalar),
    "ic_recovered_scalar": str(Integer(ic["recovered_scalar"])),
    "rho_recovered_scalar": str(Integer(rho["recovered_fixture_scalar"])),
    "generator_exact_order_verified": True,
    "ic_scalar_replay_verified": True,
    "rho_scalar_replay_verified": True,
    "same_public_point_verified": True,
    "sage_runtime_info": str(run_dir / "sage_runtime_info.json")
}
(run_dir / "independent_validation.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
print(json.dumps(result, sort_keys=True))
