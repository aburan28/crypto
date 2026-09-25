#!/usr/bin/env python3
"""Replay the frozen admission audit, then check today's narrow interface contract."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = ROOT / "research/notes/ecc2k130/rotated_solver_admission_20260925"
ORIGINAL_FILES = ("INPUT.json", "PROTOCOL.md", "audit.py", "run.py", "ci_replay.py")
NEW_FILES = (
    "PROTOCOL.md",
    "replay.py",
    "test_replay.py",
    "historical_koblitz_groebner.rs.gz",
    ".github/workflows/ecc2k130-rotated-solver-admission.yml",
)
BASE_COMMIT = "8b640f31195238242be8940c928454ad2c9f2ce6"
HISTORICAL_BUILDER = "src/cryptanalysis/koblitz_groebner.rs"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_reference(path: str) -> Path:
    rel = PurePosixPath(path)
    assert not rel.is_absolute() and ".." not in rel.parts and rel.parts
    return Path(*rel.parts)


def frozen_input() -> dict:
    own = json.loads((HERE / "FROZEN.json").read_text())
    for name in NEW_FILES:
        path = ROOT / name if name.startswith(".github/") else HERE / name
        assert sha(path.read_bytes()) == own["files"][name], name
    assert set(own["files"]) == set(NEW_FILES)
    old = json.loads((OLD / "FROZEN.json").read_text())
    assert set(old["files"]) == set(ORIGINAL_FILES)
    for name in ORIGINAL_FILES:
        assert sha((OLD / name).read_bytes()) == old["files"][name], name
    data = json.loads((OLD / "INPUT.json").read_text())
    assert data["base_commit"] == BASE_COMMIT
    assert data["domain"] == "ECC2K130-ROTATED-M56-SOLVER-ADMISSION-20260925-v1"
    assert data["reference_files"]["koblitz_builder"]["path"] == HISTORICAL_BUILDER
    return data


def historical_replay(data: dict) -> None:
    with tempfile.TemporaryDirectory(prefix="rotated-admission-historical-") as temporary:
        root = Path(temporary)
        old_copy = root / OLD.relative_to(ROOT)
        old_copy.mkdir(parents=True)
        for name in ORIGINAL_FILES + ("FROZEN.json",):
            (old_copy / name).write_bytes((OLD / name).read_bytes())
        seen = set()
        for label, item in data["reference_files"].items():
            rel = safe_reference(item["path"])
            assert rel not in seen, rel
            seen.add(rel)
            if label == "koblitz_builder":
                with gzip.open(HERE / "historical_koblitz_groebner.rs.gz", "rb") as stream:
                    raw = stream.read(1024 * 1024 + 1)
                assert len(raw) <= 1024 * 1024
            else:
                raw = (ROOT / rel).read_bytes()
            assert sha(raw) == item["sha256"], label
            target = root / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
        command = [sys.executable, str(old_copy / "ci_replay.py")]
        for suffix in ((), ("--evidence", str(OLD / "evidence/final"))):
            result = subprocess.run(command + list(suffix), cwd=root,
                                    capture_output=True, text=True, timeout=30,
                                    check=False)
            if result.returncode:
                raise AssertionError(
                    f"original replay failed ({result.returncode}): "
                    f"{result.stdout}\n{result.stderr}"
                )
            print(result.stdout.strip())
    print("HISTORICAL_REPLAY_PASS")


def function_body(source: str, name: str) -> str:
    match = re.search(rf"(?m)^pub fn {re.escape(name)}\s*\(", source)
    assert match, name
    end = source.find("\n}\n", match.end())
    assert end >= 0, name
    return source[match.start():end + 2]


def current_contract() -> dict:
    builder_path = ROOT / HISTORICAL_BUILDER
    reuse_path = ROOT / "src/cryptanalysis/polynomial_reuse.rs"
    frontend_path = ROOT / "src/cryptanalysis/koblitz_index_calculus.rs"
    builder = builder_path.read_text()
    reuse = reuse_path.read_text()
    frontend = frontend_path.read_text()
    assert re.search(r"pub const MAX_VARS:\s*usize\s*=\s*64\s*;", builder)
    body = function_body(builder, "build_decomposition_system")
    assert "basis: &[F2mElement]" in body
    assert re.search(r"if m < 2\s*\{\s*return None;", body)
    assert ".checked_mul(ell)?" in body
    assert ".checked_add((m - 2).checked_mul(n as usize)?)?" in body
    assert re.search(r"if n_vars > MAX_VARS\s*\{\s*return None;", body)
    assert "SymElement::from_subspace_vars(basis, i * ell, n, n_vars)" in body
    assert "m * ell + i * n as usize" in body
    template_start = reuse.index("impl DecompositionTemplate {")
    template_build_start = reuse.index("    pub fn build(", template_start)
    template_build_end = reuse.index("\n    /// Whether this template", template_build_start)
    template = reuse[template_build_start:template_build_end]
    assert "basis: &[F2mElement]" in template
    assert "if m < 2 || st.n == 0 || st.n > 64" in template
    assert ".checked_mul(ell)?" in template
    assert ".checked_add((m - 2).checked_mul(n as usize)?)?" in template
    assert re.search(r"if n_vars > MAX_VARS\s*\{\s*return None;", template)
    assert "SymElement::from_subspace_vars(basis, i * ell, n, n_vars)" in template
    assert "m * ell + i * n as usize" in template
    reused = function_body(reuse, "build_decomposition_system_reusing")
    assert "return build_decomposition_system(basis, x_r, b, m, st)" in reused
    assert reused.count("DecompositionTemplate::build(basis, b, m, st)") == 2
    assert "Some(template.instantiate(x_r))" in reused
    groebner = function_body(frontend, "groebner_decompose")
    assert re.search(r"let unsupported = \|\| SolveStats\s*\{\s*exhausted: true,\s*unsupported: true,", groebner)
    assert "BinaryPoint::Infinity => return (None, unsupported())" in groebner
    assert "None => return (None, unsupported())" in groebner
    assert "build_decomposition_system_reusing(" in groebner
    assert "fn pdp_admission_unsupported_is_not_a_refutation()" in frontend
    assert "fn pdp_admission_n19_m6_d2_layout_is_unsupported_not_refuted()" in frontend
    assert 5 * 2 + (5 - 2) * 13 == 49
    assert 6 * 2 + (6 - 2) * 19 == 88
    return {
        "decision": "CURRENT_GENERIC_CONTRACT_PASS",
        "scope": "generic_direct_reuse_and_groebner_frontend_only",
        "current_sha256": {
            HISTORICAL_BUILDER: sha(builder_path.read_bytes()),
            "src/cryptanalysis/polynomial_reuse.rs": sha(reuse_path.read_bytes()),
            "src/cryptanalysis/koblitz_index_calculus.rs": sha(frontend_path.read_bytes()),
        },
        "max_vars": 64,
        "shared_basis_for_all_slots": True,
        "checked_layout_arithmetic": True,
        "reuse_template_width_checked": True,
        "unsupported_is_inconclusive": True,
        "raw_bits": {"n13-m5-d2": 49, "n19-m6-d2": 88},
        "generic_builder_admits_layout": {"n13-m5-d2": True, "n19-m6-d2": False},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("freeze", "historical", "current"), required=True)
    args = parser.parse_args()
    data = frozen_input()
    if args.mode == "freeze":
        print("NEW_AND_ORIGINAL_FREEZE_PASS")
    elif args.mode == "historical":
        historical_replay(data)
    else:
        print(json.dumps(current_contract(), sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
