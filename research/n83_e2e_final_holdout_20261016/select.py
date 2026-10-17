#!/usr/bin/env python3
"""Apply the frozen N83 finalist paired selection rule to validated rows."""

import argparse
import json
import math
import statistics
from pathlib import Path

SOURCE_SHA256 = "8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12"
TARGETS_SHA256 = "407cd248797100214ff039193bb04bf89b5c2708a8cd339ca98f7aca98174062"
CURVE_PREFIX = "IC1N83Ckb1a1fb{}PDP4rootRCguidedLAgaussTDdirectISO0{}"
FINALISTS = (
    (1104, "root11", "h6f30a77abf9e"),
    (1104, "root22", "hce7cffedf9c0"),
    (1107, "root11", "h8c736ae9c7b4"),
    (1107, "root22", "hcb73613ebc5a"),
    (1112, "root11", "h49683406b22f"),
    (1112, "root22", "hbc5ec3072dc4"),
    (1120, "root11", "hc922b499a234"),
    (1120, "root22", "h82b74c23b7bc"),
)
IDS = tuple(CURVE_PREFIX.format(166 * k, suffix) for k, _, suffix in FINALISTS)
PAIRS = len(IDS) * (len(IDS) - 1) // 2
ALPHA_PER_PAIR = 0.05 / PAIRS
T_975_DF15 = 2.131449545559323
T_975_DF7 = 2.3646242515927844


def finite_positive(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value > 0


def exact_sign_flip_p(differences):
    """Two-sided paired randomization p, including equally extreme sums."""
    if len(differences) != 16 or not all(math.isfinite(x) for x in differences):
        raise ValueError("exact sign-flip test requires sixteen finite differences")
    observed = abs(math.fsum(differences))
    tolerance = 1e-12 * max(1.0, math.fsum(abs(x) for x in differences))
    extreme = 0
    for signs in range(1 << 16):
        signed = math.fsum(value if signs & (1 << i) else -value
                           for i, value in enumerate(differences))
        if abs(signed) + tolerance >= observed:
            extreme += 1
    return extreme / (1 << 16)


def paired_summary(differences):
    mean = statistics.mean(differences)
    sd = statistics.stdev(differences)
    half_width = T_975_DF15 * sd / 4
    return {"mean": mean, "median": statistics.median(differences),
            "sample_sd": sd, "lower_95": mean - half_width,
            "upper_95": mean + half_width,
            "two_sided_sign_flip_p": exact_sign_flip_p(differences)}


def admitted_matrix(data):
    if (data.get("schema") != "n83.e2e-final-validated/v1"
            or data.get("source_sha256") != SOURCE_SHA256
            or data.get("targets_sha256") != TARGETS_SHA256
            or data.get("evidence_status") != "PASS_INDEPENDENT_REPLAY"):
        return None, "SOURCE_TARGET_OR_REPLAY_GATE"
    hosts = data.get("host_strata")
    if (not isinstance(hosts, list) or len(hosts) != 2
            or {tuple(host.get("point_indices", [])) for host in hosts}
            != {tuple(range(8)), tuple(range(8, 16))}
            or any(host.get("tier_status") != "PASS_TIER1"
                   or host.get("aa_status") != "PASS_10_PERCENT_PER_PAIR"
                   for host in hosts)):
        return None, "HOST_STRATA_GATE"
    points = data.get("points")
    if (not isinstance(points, list) or len(points) != 16
            or {row.get("point_index") for row in points} != set(range(16))):
        return None, "POINT_COVERAGE_GATE"
    ordered = sorted(points, key=lambda row: row["point_index"])
    matrix = []
    for point in ordered:
        rho_by_candidate = point.get("rho_by_candidate")
        calibrations = point.get("calibration_by_candidate")
        if (not isinstance(rho_by_candidate, dict)
                or set(rho_by_candidate) != set(IDS)
                or not isinstance(calibrations, dict)
                or set(calibrations) != set(IDS)
                or any(rho.get("status") != "PASS_INDEPENDENT_REPLAY"
                       or rho.get("resource_status") != "PASS"
                       or not finite_positive(rho.get("online_ms"))
                       for rho in rho_by_candidate.values())):
            return None, "RHO_REPLAY_GATE"
        if any(item.get("status") != "PASS_10_PERCENT" for item in
               calibrations.values()):
            return None, "CALIBRATION_GATE"
        arms = point.get("arms")
        if not isinstance(arms, list) or len(arms) != len(IDS):
            return None, "IC_ARM_COVERAGE_GATE"
        by_id = {arm.get("candidate_id"): arm for arm in arms}
        if set(by_id) != set(IDS) or len(by_id) != len(arms):
            return None, "IC_CANDIDATE_ID_GATE"
        for arm in arms:
            if (arm.get("status") != "PASS_INDEPENDENT_REPLAY"
                    or arm.get("resource_status") != "PASS"
                    or not finite_positive(arm.get("one_target_online_ms"))
                    or not finite_positive(arm.get("cold_process_ms"))):
                return None, "IC_REPLAY_RESOURCE_OR_TIMING_GATE"
        matrix.append({"point_index": point["point_index"],
                       "rho_by_candidate": rho_by_candidate, "arms": by_id})
    return matrix, "PASS"


def rank_metric(matrix, metric):
    comparisons = []
    for i, first in enumerate(IDS):
        for second in IDS[i + 1:]:
            if metric == "cold_process_ms":
                differences = [point["arms"][second][metric] -
                               point["arms"][first][metric] for point in matrix]
            else:
                differences = [math.log(point["arms"][second][metric] /
                                        point["arms"][first][metric])
                               for point in matrix]
            comparisons.append({"first": first, "second": second,
                                "difference_definition": "second_minus_first" if metric == "cold_process_ms" else "log_second_over_first",
                                **paired_summary(differences)})
    winners = []
    for candidate in IDS:
        if all((comparison["mean"] > 0 if comparison["first"] == candidate
                else comparison["mean"] < 0)
               and comparison["two_sided_sign_flip_p"] <= ALPHA_PER_PAIR
               for comparison in comparisons
               if candidate in (comparison["first"], comparison["second"])):
            winners.append(candidate)
    return {"status": "SELECTED" if len(winners) == 1 else "BOUNDED_UNKNOWN",
            "selected_candidate_id": winners[0] if len(winners) == 1 else None,
            "bonferroni_alpha_per_pair": ALPHA_PER_PAIR,
            "comparisons": comparisons}


def select(data):
    matrix, admission = admitted_matrix(data)
    if matrix is None:
        return {"schema": "n83.e2e-final-selection/v1",
                "admission_status": admission,
                "cold": {"status": "BOUNDED_UNKNOWN", "selected_candidate_id": None},
                "online": {"status": "BOUNDED_UNKNOWN", "selected_candidate_id": None}}
    ratios = {candidate: [point["rho_by_candidate"][candidate]["online_ms"] /
                          point["arms"][candidate]["one_target_online_ms"]
                          for point in matrix] for candidate in IDS}
    per_host = {str(half): {candidate: values[half * 8:half * 8 + 8]
                            for candidate, values in ratios.items()}
                for half in (0, 1)}
    tier1_speedup = {}
    for half in (0, 1):
        tier1_speedup[str(half)] = {}
        for candidate in IDS:
            segment = matrix[half * 8:half * 8 + 8]
            ic = [point["arms"][candidate]["one_target_online_ms"]
                  for point in segment]
            rho = [point["rho_by_candidate"][candidate]["online_ms"]
                   for point in segment]
            differences = [reference - measured for reference, measured in zip(rho, ic)]
            mean_difference = statistics.mean(differences)
            half_width = T_975_DF7 * statistics.stdev(differences) / math.sqrt(8)
            tier1_speedup[str(half)][candidate] = {
                "status": "PASS_ALL_EIGHT_SAME_POINT" if all(x > 0 for x in differences)
                          else "NO_ALL_POINT_SPEEDUP",
                "paired_ratios": per_host[str(half)][candidate],
                "median_ratio": statistics.median(per_host[str(half)][candidate]),
                "mean_ic_online_ms": statistics.mean(ic),
                "mean_rho_online_ms": statistics.mean(rho),
                "mean_paired_rho_minus_ic_ms": mean_difference,
                "paired_difference_lower_95_ms": mean_difference - half_width,
                "paired_difference_upper_95_ms": mean_difference + half_width,
            }
    return {"schema": "n83.e2e-final-selection/v1",
            "admission_status": admission, "candidate_ids": IDS,
            "point_count": len(matrix),
            "cold": rank_metric(matrix, "cold_process_ms"),
            "online": rank_metric(matrix, "one_target_online_ms"),
            "rho_over_ic_online_by_candidate": ratios,
            "rho_over_ic_online_by_host": per_host,
            "tier1_speedup_by_host": tier1_speedup}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("validated_rows", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = select(json.loads(args.validated_rows.read_text()))
    blob = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.output:
        if args.output.exists() and args.output.read_text() != blob:
            raise ValueError("sealed selection output differs")
        args.output.write_text(blob)
    else:
        print(blob, end="")


if __name__ == "__main__":
    main()
