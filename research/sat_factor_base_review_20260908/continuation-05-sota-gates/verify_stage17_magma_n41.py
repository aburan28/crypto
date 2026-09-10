#!/usr/bin/env python3
"""Verify the source-equivalent named-generator n=41 Magma receipt."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import string


HERE = Path(__file__).resolve().parent
ARTIFACT = HERE / "stage-17-magma-n41-20260909"
CELL = "n41-l5-m3-standard-a1-f0"
SEED = 2026091301


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


STAGE16 = load_module("stage16_for_n41_magma", HERE / "verify_stage16_magma_compact.py")
STAGE15 = STAGE16.STAGE15


class VerificationError(RuntimeError):
    """The named representation or retained n=41 response is invalid."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise VerificationError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise VerificationError(f"expected JSON object in {path}")
    return value


def generator_names(variables: int) -> list[str]:
    alphabet = list(string.ascii_lowercase + string.ascii_uppercase)
    if variables <= len(alphabet):
        return alphabet[:variables]
    return alphabet + [f"z{index}" for index in range(1, variables - len(alphabet) + 1)]


def mask_indices(mask: int, variables: int) -> list[int]:
    if mask == 0:
        return []
    indices = [index for index in range(1, variables + 1) if mask & (1 << index)]
    if sum(1 << index for index in indices) != mask:
        raise VerificationError("monomial mask uses a bit outside the source-variable range")
    return indices


def render_named(variables: int, equations: list[list[int]]) -> bytes:
    names = generator_names(variables)
    polynomials = []
    for equation in equations:
        terms = []
        for mask in equation:
            indices = mask_indices(mask, variables)
            terms.append("*".join(names[index - 1] for index in indices) if indices else "1")
        polynomials.append("+".join(terms))
    prefix = (
        "SetGPU(false);SetSeed(1);rr<"
        + ",".join(names)
        + f'>:=BooleanPolynomialRing({variables},"grevlex");ff:=['
        + ",".join(polynomials)
        + "];"
    )
    suffix = (
        'ii:=ideal<rr|ff>;cc:=Cputime();ww:=Realtime();gg,dd:=GroebnerBasis(ii:Al:="Direct",'
        'Faugere:=true,Dense:=false,Nthreads:=1);cc:=Cputime(cc);ww:=Realtime(ww);printf '
        '"KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1\\n";printf '
        '"KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse\\n";if #gg eq 1 and gg[1] eq rr!1 '
        'then printf "KOBLITZ_MAGMA_STATUS=UNSAT\\n";else printf '
        '"KOBLITZ_MAGMA_STATUS=SAT\\n";end if;printf '
        '"KOBLITZ_MAGMA_F4_DEGREES=%o\\n",dd;printf '
        '"KOBLITZ_MAGMA_BASIS_SIZE=%o\\n",#gg;printf '
        '"KOBLITZ_MAGMA_CPU_SECONDS=%o\\n",cc;printf '
        '"KOBLITZ_MAGMA_WALL_SECONDS=%o\\n",ww;quit;'
    )
    return (prefix + suffix).encode()


def named_inventory(equivalence: dict) -> dict:
    rows = []
    for row in equivalence["rows"]:
        anf_path, _ = STAGE16.task_paths(row["seed"], row["cell_id"])
        variables, equations = STAGE16.parse_anf(anf_path)
        rendered = render_named(variables, equations)
        rows.append(
            {
                "seed": row["seed"],
                "cell_id": row["cell_id"],
                "bytes": len(rendered),
                "sha256": sha256_bytes(rendered),
                "fits_calculator": len(rendered) <= 50_000,
            }
        )
    eligible = [
        {"seed": row["seed"], "cell_id": row["cell_id"]}
        for row in rows
        if row["fits_calculator"]
    ]
    if len(rows) != 20 or len(eligible) != 15:
        raise VerificationError("named-generator eligibility inventory changed")
    if any(item["cell_id"].startswith("n59-") for item in eligible):
        raise VerificationError("an n=59 input unexpectedly fits the calculator")
    return {"tasks": len(rows), "eligible_tasks": eligible, "rows": rows}


def summarize() -> dict:
    expected_names = {
        "seed-2026091301-named-input.magma",
        "seed-2026091301-equivalence.json",
        "seed-2026091301-response.xml",
        "stage-17-summary.json",
    }
    artifacts = list(ARTIFACT.iterdir())
    if {path.name for path in artifacts} != expected_names or any(
        path.is_symlink() or not path.is_file() for path in artifacts
    ):
        raise VerificationError("Stage 17 artifact inventory changed")
    equivalence_inventory = STAGE16.complete_equivalence_inventory()
    inventory = named_inventory(equivalence_inventory)
    source = next(
        row
        for row in equivalence_inventory["rows"]
        if row["seed"] == SEED and row["cell_id"] == CELL
    )
    anf_path, _ = STAGE16.task_paths(SEED, CELL)
    variables, equations = STAGE16.parse_anf(anf_path)
    named = render_named(variables, equations)
    input_path = ARTIFACT / "seed-2026091301-named-input.magma"
    if input_path.read_bytes() != named:
        raise VerificationError("retained n=41 named input differs from deterministic rendering")
    equivalence = read_json(ARTIFACT / "seed-2026091301-equivalence.json")
    expected_equivalence = {
        "seed": SEED,
        "cell_id": CELL,
        "source_anf_sha256": source["anf_sha256"],
        "source_magma_sha256": source["verbose_magma_sha256"],
        "named_input_sha256": sha256_bytes(named),
        "named_input_bytes": len(named),
        "equations": source["equations"],
        "variables": source["variables"],
        "equation_monomial_lists_equal": True,
        "representation": "one named Boolean-ring generator per source variable",
    }
    if equivalence != expected_equivalence:
        raise VerificationError("n=41 equivalence receipt changed")
    response = STAGE15.parse_calculator_xml(ARTIFACT / "seed-2026091301-response.xml")
    STAGE15.require_clean_f4_output(response["output"], SEED)
    if (
        response["service"]["warning"] is not None
        or response["service"]["alert"] is not None
    ):
        raise VerificationError("n=41 response contains a service diagnostic")
    terminal = STAGE15.parse_magma_terminal(response["output"])
    if terminal is None or terminal["terminal_status"] != "sat":
        raise VerificationError("n=41 response lacks a proper-ideal terminal")
    if terminal["f4_step_degrees"] != [2, 3, 3, 3, 3] or terminal["basis_size"] != 43:
        raise VerificationError("n=41 F4 terminal identity changed")
    return {
        "schema": "koblitz_magma_n41_stage17_summary.v1",
        "evidence_class": "source_equivalent_external_service_receipt",
        "source_equivalence": equivalence,
        "named_representation_inventory": inventory,
        "terminal": terminal,
        "service": response["service"],
        "response_sha256": response["response_sha256"],
        "classification": "sat_basis_certificate_unverified_model",
        "validated_magma_point_witnesses": 0,
        "magma_process_resource_receipts_complete": False,
        "full_magma_matrix_executed": False,
        "full_solver_matrix_gate_passed": False,
        "claim": (
            "One source-equivalent n=41 standard instance has a calculator-reported direct-F4 "
            "proper-ideal terminal. It lacks a point witness and process-scoped resource receipt; "
            "this is partial PDP evidence, not the repeated Magma matrix or a SOTA result."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert generator_names(3) == ["a", "b", "c"]
        assert mask_indices(6, 3) == [1, 2]
        print(json.dumps({"self_test": "pass"}, indent=2))
        return
    summary = summarize()
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.expected is not None and args.expected.read_text() != rendered:
        raise VerificationError("recomputed Stage 17 summary differs from expected summary")
    if args.output is not None:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
