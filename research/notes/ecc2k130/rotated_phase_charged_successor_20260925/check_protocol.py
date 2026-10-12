#!/usr/bin/env python3
"""Hash-only and metadata-only validation of the charged successor design.

This checker deliberately does not open phase outcomes or run a PDP solver.
The executable implementation requires its own later exact-head freeze.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_factors(raw: bytes, slots: int, choices: int, label: str) -> None:
    factors = json.loads(raw)
    require(len(factors) == slots, f"{label} slot count drift")
    for slot in factors:
        require(len(slot) == choices, f"{label} factor size drift")
        require(all(isinstance(point, list) and len(point) == 2 and
                    all(isinstance(coord, int) for coord in point)
                    for point in slot), f"{label} nonfinite factor choice")
        require(slot.count([0, 1]) == 1, f"{label} four-torsion factor drift")


def main() -> None:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    require(frozen["status"] == "design_freeze_no_solver_outcome", "freeze status drift")
    for relative, expected in frozen["source_sha256"].items():
        require(digest((ROOT / relative).read_bytes()) == expected, f"{relative} source SHA drift")
    config = json.loads((HERE / "INPUT.json").read_text())
    require(config["status"] == "design_freeze_no_solver_outcome", "bad status")
    require(config["schema"] == "ecc2k130_rotated_phase_charged_successor_design_v1", "bad schema")
    require((config["n"], config["m"], config["d"], config["q"]) == (19, 6, 2, 130873), "bad field or group")
    require(config["phase_order_k8"] == [(5 * i) % 19 for i in range(8)], "phase order drift")
    require(config["phase_prefix_k4"] == config["phase_order_k8"][:4], "prefix drift")
    require(config["beta_order_k5"] == [3, 338435, 303097, 464276, 42605], "beta order drift")
    require(config["torsion_order"] == [[1, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1]], "torsion drift")
    require(config["software_gate"]["required_parent_prs_merged_before_run"] == [802, 804, 807], "parent gate drift")
    require(config["software_gate"]["dimacs_id_max"] == 2**31 - 1, "DIMACS width drift")
    require(config["resource_gates"]["query_cnf_byte_max"] == 256 * 2**20, "CNF cap drift")
    require(config["resource_gates"]["proof_byte_max"] == 2**30, "proof cap drift")
    require(config["rho_seed_repetitions"] == [0, 1, 2], "rho seeds drift")
    require(config["single_beta_control_phase"] == 0, "single beta control drift")
    require(config["factor_identity_choices_per_slot"] == 0, "factor O-choice drift")

    selected: dict[str, set[tuple[int, int]]] = {}
    for key, prefix, count in [("training", "tr", 256), ("holdout", "ho", 64)]:
        item = config[key]
        path = ROOT / item["path"]
        raw = path.read_bytes()
        require(digest(raw) == item["sha256"], f"{key} SHA drift")
        doc = json.loads(raw)
        cases = doc["targets"]
        require(len(cases) == count == item["case_count"], f"{key} count drift")
        require([r["case_id"] for r in cases] == [f"{prefix}-{i:03d}" for i in range(count)], f"{key} order drift")
        if key == "holdout":
            require(all(set(r) == {"case_id", "point"} for r in cases), "holdout labels leaked")
        selected[key] = {tuple(r["point"]) for r in cases}
        require(len(selected[key]) == count, f"duplicate {key} point")
    require(not selected["training"] & selected["holdout"], "training/holdout overlap")

    refs = config["source_archives"]
    archive = ROOT / "research/notes/ecc2k130/rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
    require(digest(archive.read_bytes()) == refs["n13_n19_corpus_sha256"], "corpus SHA drift")
    with tarfile.open(archive, "r:gz") as tf:
        for member, key in [
            ("raw/n13-m5/targets.json", "n13_targets_sha256"),
            ("raw/n19-m6/targets.json", "n19_targets_sha256"),
            ("raw/n13-m5/factors.json", "n13_factors_sha256"),
            ("raw/n19-m6/factors.json", "n19_factors_sha256"),
        ]:
            source = tf.extractfile(member)
            require(source is not None, f"{member} absent")
            raw = source.read()
            require(digest(raw) == refs[key], f"{member} SHA drift")
            if member.endswith("factors.json"):
                check_factors(raw, 5 if "n13" in member else 6,
                              5 if "n13" in member else 7, member)
    portfolio = ROOT / "research/notes/ecc2k130/rotated_beta_sweep_20260925/evidence/raw.tar.gz"
    require(digest(portfolio.read_bytes()) == refs["additional_beta_archive_sha256"], "other-beta SHA drift")
    with tarfile.open(portfolio, "r:gz") as tf:
        for beta in config["beta_order_k5"][1:]:
            member = f"raw/beta-{beta}/factors.json"
            source = tf.extractfile(member)
            require(source is not None, f"{member} absent")
            check_factors(source.read(), 6, 7, member)

    require((6 * (2 * 19 + 1) + 4 * (2 * 19 + 1) + (2 * 19 + 1) + 5 * 19) == 524, "layout drift")
    print("PASS: frozen design inputs, target split, source archives and width/cap literals")


if __name__ == "__main__":
    main()
