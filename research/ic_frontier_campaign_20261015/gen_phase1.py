"""Write the Phase 1 specs from the frozen per-size parameters (params.json).

Lays out specs only; computes no result.  One spec per curve, every arm of the
family on every curve, hold-out public targets (seed 20261015200), 16 targets
x 2 rounds, 1 warm-up, counts only.

    python3 gen_phase1.py        # writes phase1/spec-*.json
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HERE, "params.json")))
OUT = os.path.join(HERE, "phase1")
os.makedirs(OUT, exist_ok=True)

MEAS = {"rounds": 2, "warmup": 1, "order": "alternate", "seed": 20261015201,
        "isolation_required": "L0", "timeout_seconds": 3600}
WL = {"targets_per_curve": 16, "target_seed": 20261015200, "target_kind": "public"}


def ic(name, fb, oracle, trials, role="candidate"):
    return {"name": name, "role": role, "method": {"id": "ic.pipeline", "params": {
        "factor_base": fb, "oracle": oracle, "max_trials": str(trials)}}}


for k in P["koblitz"]:
    arms = [
        {"name": "rho-strong", "role": "reference", "method": {"id": "rho.signed_frobenius_strong"}},
        {"name": "rho-frob", "role": "baseline", "method": {"id": "rho.signed_frobenius"}},
        {"name": "bsgs-neg", "role": "baseline", "method": {"id": "bsgs.negation"}},
        ic("ic-compact", f"compact-orbit-scan:columns={k['compact_columns']},raw_x_cap=1000000",
           "mitm-frobenius-counted:m=2", 50000000),
        ic("ic-negfold", f"koblitz-standard-subspace:dimension={k['negfold_dimension']}",
           "pdp3-koblitz:m=2,engine=rootfind-negfold", 50000000),
    ]
    if k["n"] <= 23:
        arms.append(ic("ic-f4", "koblitz-standard-subspace:dimension=5",
                       "pdp3-koblitz:m=3,engine=inherited-f4,degree=3,node_budget=8192", 200000,
                       role="baseline"))
    arms.append({"name": "ic-shared", "role": "candidate", "method": {"id": "ic.shared_rank", "params": {
        "factor_base": "compact-orbit-scan:columns=8,raw_x_cap=1000000",
        "rank_seed": "20261015202", "rank_max_trials": "1000000", "target_max_attempts": "512"}}})
    arms.append({"name": "rho-strong-aa", "role": "control", "method": {"id": "rho.signed_frobenius_strong"}})
    spec = {"schema": "ecbench.spec/v1",
            "label": f"IC frontier campaign Phase 1: Koblitz a={k['a']} n={k['n']}",
            "description": "Phase 1 counted frontier (research/ic_frontier_campaign_20261015/PROTOCOL.md, "
                           "frozen): every Koblitz arm on hold-out public targets, parameters from "
                           "params.json chosen on the train panel. Counts only; ranking by the AGENTS.md "
                           "5a band. Public synthetic; no security claim.",
            "workloads": dict(WL, curves=[{"kind": "koblitz", "a": k["a"], "n": k["n"]}]),
            "arms": arms, "measurement": MEAS}
    json.dump(spec, open(os.path.join(OUT, f"spec-K{k['n']}.json"), "w"), indent=1)

curves = [json.loads(l) for l in open(os.path.join(HERE, "prime_curves.jsonl"))]
for c, p in zip(curves, P["prime"]):
    r = c["r"]
    N = p["pair_size"]
    arms = [
        {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
        {"name": "bsgs-neg", "role": "baseline", "method": {"id": "bsgs.negation"}},
        {"name": "kangaroo", "role": "baseline", "method": {"id": "kangaroo.vow"}},
        ic("ic-pair", f"prime-abscissa:size={N}", "mitm:negation_folded=1", int(40 * r / N) + 100000),
        {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}},
    ]
    spec = {"schema": "ecbench.spec/v1",
            "label": f"IC frontier campaign Phase 1: {p['bits']}-bit a=-3 prime-order curve",
            "description": "Phase 1 counted frontier (PROTOCOL.md, frozen): generic arms and the two-summand "
                           "pair-table IC with the train-panel base size. Counts only; AGENTS.md 5a band. "
                           "Public synthetic; no security claim.",
            "workloads": dict(WL, curves=[c]), "arms": arms, "measurement": MEAS}
    json.dump(spec, open(os.path.join(OUT, f"spec-P{p['bits']}.json"), "w"), indent=1)

pkm = [json.loads(l) for l in open(os.path.join(HERE, "..", "pkm_tower_ecbench_20261013", "curves.jsonl"))]
for c, (bits, t) in zip(pkm, P["pkm"]):
    arms = [
        {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
        ic("ic-pkm", f"prime-kummer:t={t}", "rootfind-mitm:m=3", 2000000),
        {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}},
    ]
    spec = {"schema": "ecbench.spec/v1",
            "label": f"IC frontier campaign Phase 1: PKM tower ladder, {bits}-bit",
            "description": "Phase 1 (PROTOCOL.md): the Kummer-tower arm on its registered 2-adic ladder "
                           "(research/pkm_tower_ecbench_20261013), t as registered there. Counts only.",
            "workloads": dict(WL, curves=[c]), "arms": arms, "measurement": MEAS}
    json.dump(spec, open(os.path.join(OUT, f"spec-PKM{bits}.json"), "w"), indent=1)
print("ok")
