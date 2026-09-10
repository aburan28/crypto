#!/usr/bin/env python3
"""Verify compact source-equivalent Magma F4 calculator receipts."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
PANEL = HERE / "stage-13-panel-20260909"
ARTIFACT = HERE / "stage-16-magma-compact-20260909"
PROTOCOL = HERE / "stage-13-pdp-panel-protocol.json"
CELL = "n31-l5-m3-ggmp-a0-f0"
SEEDS = (2026091303, 2026091305)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


STAGE15 = load_module("stage15_for_compact_magma", HERE / "verify_stage15_magma_calculator.py")
STAGE13 = load_module("stage13_for_compact_magma", HERE / "verify_stage13_pdp_panel.py")


class VerificationError(RuntimeError):
    """A compact representation or retained response violates its binding."""


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


def parse_anf(path: Path) -> tuple[int, list[list[int]]]:
    lines = path.read_text().splitlines()
    header = lines[0].split() if lines else []
    if len(header) != 4 or header[:2] != ["p", "cnf"]:
        raise VerificationError(f"{path}: malformed ANF header")
    try:
        variables, expected_equations = map(int, header[2:])
    except ValueError as error:
        raise VerificationError(f"{path}: malformed ANF dimensions") from error
    equations: list[list[int]] = []
    for line_number, line in enumerate(lines[1:], 2):
        tokens = line.split()
        if len(tokens) < 3 or tokens[0] != "x" or tokens[-1] != "0":
            raise VerificationError(f"{path}:{line_number}: malformed XOR row")
        index = 1
        masks: list[int] = []
        has_true_token = False
        while tokens[index] != "0":
            token = tokens[index]
            if token == "T":
                if has_true_token:
                    raise VerificationError(f"{path}:{line_number}: duplicate T token")
                has_true_token = True
                index += 1
                continue
            if token.startswith("."):
                try:
                    degree = int(token[1:])
                except ValueError as error:
                    raise VerificationError(f"{path}:{line_number}: bad degree token") from error
                names = tokens[index + 1 : index + 1 + degree]
                if degree < 2 or len(names) != degree:
                    raise VerificationError(f"{path}:{line_number}: truncated monomial")
                try:
                    indices = [int(name) for name in names]
                except ValueError as error:
                    raise VerificationError(f"{path}:{line_number}: bad variable") from error
                index += degree + 1
            else:
                try:
                    indices = [int(token)]
                except ValueError as error:
                    raise VerificationError(f"{path}:{line_number}: bad variable") from error
                index += 1
            if (
                len(set(indices)) != len(indices)
                or any(variable < 1 or variable > variables for variable in indices)
            ):
                raise VerificationError(f"{path}:{line_number}: invalid Boolean monomial")
            mask = 0
            for variable in indices:
                mask |= 1 << variable
            masks.append(mask)
        # WDSat XOR rows have an implicit true right-hand side.  Its `T`
        # token toggles that constant.  The exact comparison with the frozen
        # Magma list below independently checks this conversion for all rows.
        if not has_true_token:
            masks.append(0)
        if len(set(masks)) != len(masks):
            raise VerificationError(f"{path}:{line_number}: duplicate Boolean monomial")
        equations.append(masks)
    if len(equations) != expected_equations:
        raise VerificationError(f"{path}: ANF equation count changed")
    return variables, equations


def parse_verbose_magma(path: Path, variables: int) -> list[list[int]]:
    text = path.read_text()
    ring = f'R := BooleanPolynomialRing({variables}, "grevlex");\n'
    if ring not in text:
        raise VerificationError(f"{path}: Boolean ring declaration changed")
    try:
        body = text.split("F := [", 1)[1].split("];\nI := ideal<R | F>;", 1)[0]
    except IndexError as error:
        raise VerificationError(f"{path}: cannot isolate Magma F list") from error
    equations = []
    for polynomial in body.split(",\n  "):
        masks = []
        for term in polynomial.strip().split(" + "):
            if term == "1":
                masks.append(0)
                continue
            indices = [int(value) for value in re.findall(r"X\[(\d+)\]", term)]
            if not indices or term != "*".join(f"X[{value}]" for value in indices):
                raise VerificationError(f"{path}: unsupported Magma term {term!r}")
            if len(set(indices)) != len(indices) or any(
                variable < 1 or variable > variables for variable in indices
            ):
                raise VerificationError(f"{path}: invalid Magma monomial")
            mask = 0
            for variable in indices:
                mask |= 1 << variable
            masks.append(mask)
        if len(set(masks)) != len(masks):
            raise VerificationError(f"{path}: duplicate Magma monomial")
        equations.append(masks)
    return equations


def format_mask(mask: int) -> str:
    return min(str(mask), hex(mask), key=len)


def render_compact(variables: int, equations: list[list[int]]) -> bytes:
    prefix = (
        f'SetGPU(false);SetSeed(1);R:=BooleanPolynomialRing({variables},"grevlex");'
        "P:=func<Q|BooleanPolynomial(R,Q)>;F:=["
        + ",".join(
            "P([" + ",".join(format_mask(mask) for mask in equation) + "])"
            for equation in equations
        )
        + "];"
    )
    suffix = (
        'I:=ideal<R|F>;c:=Cputime();w:=Realtime();G,D:=GroebnerBasis(I:Al:="Direct",'
        'Faugere:=true,Dense:=false,Nthreads:=1);c:=Cputime(c);w:=Realtime(w);printf '
        '"KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1\\n";printf '
        '"KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse\\n";if #G eq 1 and G[1] eq R!1 '
        'then printf "KOBLITZ_MAGMA_STATUS=UNSAT\\n";else printf '
        '"KOBLITZ_MAGMA_STATUS=SAT\\n";end if;printf '
        '"KOBLITZ_MAGMA_F4_DEGREES=%o\\n",D;printf '
        '"KOBLITZ_MAGMA_BASIS_SIZE=%o\\n",#G;printf '
        '"KOBLITZ_MAGMA_CPU_SECONDS=%o\\n",c;printf '
        '"KOBLITZ_MAGMA_WALL_SECONDS=%o\\n",w;quit;'
    )
    return (prefix + suffix).encode()


def task_paths(seed: int, cell_id: str) -> tuple[Path, Path]:
    instance = PANEL / "tasks" / f"seed-{seed}" / cell_id / "matrix" / cell_id
    return instance / "instance.anf", instance / "instance.magma"


def verify_constructor_sanity() -> dict:
    expected_input = '''R := BooleanPolynomialRing(3, "grevlex");
X := [R.i : i in [1..3]];
verbose := X[1]*X[2] + X[3] + 1;
compact := BooleanPolynomial(R, [6, 8, 0]);
printf "KOBLITZ_MAGMA_MASK_SCHEMA=koblitz_magma_boolean_mask_sanity.v1\\n";
printf "KOBLITZ_MAGMA_MASK_EQUAL=%o\\n", verbose eq compact;
printf "KOBLITZ_MAGMA_MASK_VERBOSE=%o\\n", verbose;
printf "KOBLITZ_MAGMA_MASK_COMPACT=%o\\n", compact;
quit;
'''.encode()
    input_path = ARTIFACT / "boolean-mask-sanity.magma"
    if input_path.read_bytes() != expected_input:
        raise VerificationError("Boolean mask sanity input changed")
    response = STAGE15.parse_calculator_xml(ARTIFACT / "boolean-mask-sanity.xml")
    expected_output = (
        "KOBLITZ_MAGMA_MASK_SCHEMA=koblitz_magma_boolean_mask_sanity.v1\n"
        "KOBLITZ_MAGMA_MASK_EQUAL=true\n"
        "KOBLITZ_MAGMA_MASK_VERBOSE=$.1*$.2 + $.3 + 1\n"
        "KOBLITZ_MAGMA_MASK_COMPACT=$.1*$.2 + $.3 + 1\n\n"
    )
    if (
        response["service"]["warning"] is not None
        or response["service"]["alert"] is not None
        or response["output"] != expected_output
    ):
        raise VerificationError("Boolean mask constructor sanity response changed")
    return {
        "input_path": input_path.name,
        "input_sha256": sha256_bytes(expected_input),
        "response_path": "boolean-mask-sanity.xml",
        "response_sha256": response["response_sha256"],
        "verbose_equals_compact": True,
        "variable_mask_convention": "variable i is represented by integer bit 1 << i",
    }


def complete_equivalence_inventory() -> dict:
    protocol = read_json(PROTOCOL)
    STAGE13.validate_protocol(protocol)
    panel_summary = STAGE13.summarize(
        protocol,
        PANEL,
        write_receipts=False,
        allow_incomplete=False,
    )
    if (
        panel_summary.get("verified_tasks") != 20
        or panel_summary.get("panel_artifact_complete") is not True
    ):
        raise VerificationError("Stage 13 archive is not the custody-verified frozen panel")
    seeds = protocol.get("replicate_seeds")
    cells = protocol.get("cells")
    if not isinstance(seeds, list) or not isinstance(cells, list):
        raise VerificationError("frozen protocol lacks seeds or cells")
    rows = []
    for seed in seeds:
        for cell in cells:
            cell_id = cell["id"]
            anf_path, magma_path = task_paths(seed, cell_id)
            variables, anf = parse_anf(anf_path)
            verbose = parse_verbose_magma(magma_path, variables)
            if anf != verbose:
                raise VerificationError(f"seed {seed}/{cell_id}: ANF and Magma monomials differ")
            binding = read_json(
                PANEL / "tasks" / f"seed-{seed}" / cell_id / "source-binding.json"
            )
            for name, path in (
                ("wdsat_anf", anf_path),
                ("magma_boolean_f4", magma_path),
            ):
                frozen = binding.get("exports", {}).get(name)
                data = path.read_bytes()
                if not isinstance(frozen, dict) or (
                    frozen.get("bytes") != len(data)
                    or frozen.get("sha256") != sha256_bytes(data)
                ):
                    raise VerificationError(
                        f"seed {seed}/{cell_id}: {name} differs from source binding"
                    )
            compact = render_compact(variables, anf)
            rows.append(
                {
                    "seed": seed,
                    "cell_id": cell_id,
                    "basis": cell["basis"],
                    "variables": variables,
                    "equations": len(anf),
                    "anf_sha256": sha256_bytes(anf_path.read_bytes()),
                    "verbose_magma_sha256": sha256_bytes(magma_path.read_bytes()),
                    "compact_bytes": len(compact),
                    "compact_sha256": sha256_bytes(compact),
                    "fits_calculator": len(compact) <= 50_000,
                }
            )
    if len(rows) != 20:
        raise VerificationError("compact equivalence inventory is not the frozen 20 tasks")
    return {
        "tasks": len(rows),
        "all_equation_monomial_lists_equal": True,
        "calculator_eligible": [
            {"seed": row["seed"], "cell_id": row["cell_id"]}
            for row in rows
            if row["fits_calculator"]
        ],
        "rows": rows,
    }


def verify_case(seed: int, inventory: dict) -> dict:
    row = next(
        item
        for item in inventory["rows"]
        if item["seed"] == seed and item["cell_id"] == CELL
    )
    stored_input = ARTIFACT / f"seed-{seed}-compact-input.magma"
    compact = stored_input.read_bytes()
    if compact != render_compact(*parse_anf(task_paths(seed, CELL)[0])):
        raise VerificationError(f"seed {seed}: retained compact input changed")
    equivalence = read_json(ARTIFACT / f"seed-{seed}-equivalence.json")
    expected_equivalence = {
        "seed": seed,
        "cell_id": CELL,
        "source_anf_sha256": row["anf_sha256"],
        "source_magma_sha256": row["verbose_magma_sha256"],
        "compact_input_sha256": row["compact_sha256"],
        "compact_input_bytes": row["compact_bytes"],
        "equations": row["equations"],
        "variables": row["variables"],
        "equation_monomial_lists_equal": True,
        "constant_mapping": (
            "WDSat XOR rows without T have RHS true and map to a Magma constant term"
        ),
    }
    if equivalence != expected_equivalence:
        raise VerificationError(f"seed {seed}: equivalence receipt changed")
    response = STAGE15.parse_calculator_xml(ARTIFACT / f"seed-{seed}-response.xml")
    STAGE15.require_clean_f4_output(response["output"], seed)
    if (
        response["service"]["warning"] is not None
        or response["service"]["alert"] is not None
    ):
        raise VerificationError(f"seed {seed}: calculator response contains a warning")
    terminal = STAGE15.parse_magma_terminal(response["output"])
    if terminal is None or terminal["terminal_status"] != "sat":
        raise VerificationError(f"seed {seed}: missing proper-ideal terminal")
    if terminal["f4_step_degrees"] != [2, 2, 3, 3, 4, 3, 4] or terminal["basis_size"] != 87:
        raise VerificationError(f"seed {seed}: unexpected GGMP F4 terminal identity")
    return {
        "seed": seed,
        "cell_id": CELL,
        "representation": "compact BooleanPolynomial monomial-mask list",
        "source_equivalence": equivalence,
        "input_path": stored_input.name,
        "response_path": f"seed-{seed}-response.xml",
        "terminal": terminal,
        "service": response["service"],
        "response_sha256": response["response_sha256"],
        "classification": "sat_basis_certificate_unverified_model",
        "process_scoped_resources_complete": False,
        "point_witness_validated": False,
    }


def summarize() -> dict:
    expected_names = {
        "boolean-mask-sanity.magma",
        "boolean-mask-sanity.xml",
        "stage-16-summary.json",
        *(f"seed-{seed}-{suffix}" for seed in SEEDS for suffix in (
            "compact-input.magma",
            "equivalence.json",
            "response.xml",
        )),
    }
    artifacts = list(ARTIFACT.iterdir())
    if {path.name for path in artifacts} != expected_names or any(
        path.is_symlink() or not path.is_file() for path in artifacts
    ):
        raise VerificationError("Stage 16 artifact inventory changed")
    expected_responses = {
        ARTIFACT / "boolean-mask-sanity.xml",
        *(ARTIFACT / f"seed-{seed}-response.xml" for seed in SEEDS),
    }
    if set(ARTIFACT.glob("*.xml")) != expected_responses:
        raise VerificationError("retained compact-response inventory changed")
    constructor_sanity = verify_constructor_sanity()
    inventory = complete_equivalence_inventory()
    cases = [verify_case(seed, inventory) for seed in SEEDS]
    if STAGE15.terminal_identity(cases[0]["terminal"]) != STAGE15.terminal_identity(
        cases[1]["terminal"]
    ):
        raise VerificationError("compact GGMP terminal identities disagree")
    eligible_ggmp = [
        row
        for row in inventory["rows"]
        if row["basis"] == "ggmp" and row["fits_calculator"]
    ]
    if [(row["seed"], row["cell_id"]) for row in eligible_ggmp] != [
        (seed, CELL) for seed in SEEDS
    ]:
        raise VerificationError("additional compact GGMP eligibility set changed")
    return {
        "schema": "koblitz_magma_compact_stage16_summary.v1",
        "evidence_class": "source_equivalent_external_service_receipt",
        "constructor_sanity": constructor_sanity,
        "equivalence_inventory": inventory,
        "additional_compact_ggmp_responses": len(cases),
        "validated_magma_point_witnesses": 0,
        "magma_process_resource_receipts_complete": False,
        "full_magma_matrix_executed": False,
        "full_solver_matrix_gate_passed": False,
        "cases": cases,
        "claim": (
            "Two additional n=31 GGMP systems have source-equivalent compact direct-F4 "
            "proper-ideal terminals. The compact renderer is monomial-list equivalent to all "
            "twenty frozen ANF/Magma pairs. These receipts lack point witnesses and process-"
            "scoped resources and do not complete the Magma matrix or support a SOTA claim."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert format_mask(8) == "8"
        assert format_mask(65536) == "65536"
        assert format_mask(1 << 60) == "0x1000000000000000"
        print(json.dumps({"self_test": "pass"}, indent=2))
        return
    summary = summarize()
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.expected is not None and args.expected.read_text() != rendered:
        raise VerificationError("recomputed Stage 16 summary differs from expected summary")
    if args.output is not None:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
