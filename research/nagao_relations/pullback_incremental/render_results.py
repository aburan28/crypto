"""Render a static result note from a frozen comparison, without rerunning solvers."""
import argparse
import json
from pathlib import Path


def render(data):
    text = ['# Incremental pullback benchmark results', '',
        'Engineering candidate; calibrated normalized-cost classification remains pending. '
        'No end-to-end speedup is established. Counts below come from [comparison.json](comparison.json); '
        'the derivation, frozen gate and limitations are in [README.md](README.md).', '',
        f"The run contains {data['stage_run_count']} solver cells and {data['dlp_run_count']} verified cold ECDLP runs. "
        f"The audit checked {data['dlp_oracle_checks']} decomposition calls across "
        f"{data['dlp_distinct_oracle_targets']} distinct targets, including empty answers. "
        f"It reproduced {data['stage_baseline_exact_counter_replays']} completed prior stage-counter records and "
        f"all {data['dlp_baseline_exact_counter_replays']} prior cold baseline-counter records exactly. "
        f"{len(data['stage_baseline_censored_replays'])} formerly completed stage replay(s) were censored by this run's deadline; "
        'they are retained and excluded from equal-output ratios.', '',
        '| Variant / reference | Full-DLP S | Cost/rho | Cost/floor | Class/status |',
        '|---|---|---|---|---|']
    for variant in ('coefficient-pullback','normalized-pullback','cached-pullback','s4-symmetric','chained-s3','rho'):
        status = 'Engineering candidate; pending calibration' if variant in ('normalized-pullback','cached-pullback') else 'Reference/control'
        text.append(f'| {variant} | null | null | null | {status} |')
    text += ['', '## Complete eleven-bit enumeration gate', '',
        'Each corpus has four targets. These pairs complete with identical outputs. '
        'The coordinate cache preserves the normalized field vector exactly. The gate concerns '
        'multiplications; the API sum is an uncalibrated diagnostic.', '',
        '| Corpus | Variant | Additions | Multiplications | Squarings | Field-API sum |',
        '|---|---|---:|---:|---:|---:|']
    for gate in data['stage_diagnostic_gates']:
        for variant, counts in gate['field_api_counts'].items():
            text.append(f"| {gate['corpus']} | {variant} | {counts['additions']:,} | {counts['multiplications']:,} | {counts['squarings']:,} | {counts['fieldOperations']:,} |")
    for gate in data['stage_diagnostic_gates']:
        text += ['', f"{gate['corpus'].capitalize()}: {100*gate['multiplication_reduction_fraction']:.2f}% fewer multiplications; "
                 f"{gate['complete_equal_output_targets']}/4 complete equal-output targets. Frozen 10% diagnostic gate: "
                 f"{'PASS' if gate['diagnostic_gate_passed'] else 'FAIL'}."]
    text += ['', '## Resolution at the fixed budget', '',
        'Each row includes four targets in first mode plus the same four in enumeration mode. '
        'Resolved includes proved-empty cases. Relations found during incomplete enumeration do '
        'not make that cell resolved. Supported targets and uniform targets are distinct strata '
        'in the raw file; their mixture is not a natural relation-yield estimate.', '',
        '| Corpus | n / d | Variant | Resolved / 8 | Resolved within 3 s / 8 |',
        '|---|---|---|---:|---:|']
    keys = sorted({(r['corpus'],r['n'],r['d'],r['variant']) for r in data['stage_groups']})
    for corpus,n,d,variant in keys:
        rr=[r for r in data['stage_groups'] if (r['corpus'],r['n'],r['d'],r['variant'])==(corpus,n,d,variant)]
        text.append(f"| {corpus} | {n} / {d} | {variant} | {sum(r['resolved'] for r in rr)} | {sum(r['resolved_within_budget'] for r in rr)} |")
    text += ['', '## Cold ECDLP field-operation vectors', '',
        'Each row sums three cold runs on the indicated corpus and group. Every run rebuilds '
        'setup and verifies its scalar. Scalar modular counts and phase times are separate in '
        'the comparison file. These partial vectors are not calibrated total operations. '
        'Times are descriptive; one repetition cannot establish a paired runtime-confidence claim.', '',
        '| Corpus | Field bits / subgroup | Variant | Additions | Multiplications | Squarings | Cold seconds |',
        '|---|---|---|---:|---:|---:|---:|']
    for row in data['dlp_groups']:
        c=row['field_api_counts']
        text.append(f"| {row['corpus']} | {row['n']} / {row['subgroup_order']} | {row['variant']} | {c['additions']:,} | {c['multiplications']:,} | {c['squarings']:,} | {row['cold_seconds']:.6f} |")
    equal=sum(row['identical_attempt_path'] for row in data['dlp_attempt_paths'])
    text += ['', f"{equal}/{len(data['dlp_attempt_paths'])} candidate/reference pairs followed identical full attempt paths. "
        'All costs remain in the totals if a valid first relation changes the downstream path. '
        'Fresh seeds can repeat Q in tiny groups: five-bit fresh seeds provide one new distinct '
        'target beyond the frozen panel; nine-bit fresh seeds provide two. The raw runs and '
        'distinctness audit preserve these collisions.', '',
        'Rho is the same-group, same-target control. No field/scalar/control conversion or '
        'calibrated full-cost floor has been measured, so S and boundary ratios remain null. '
        'The branch exponent and support-counting ceiling are unchanged.', '',
        '## Evidence and next gate', '',
        'The initial harness import failure is preserved in `results/`; it occurred before '
        'any solver measurement. The accepted campaign is `results_v2/`. Its raw records '
        'include targets, partial outputs, failures, phase costs, scalar certificates and source '
        'hashes. Neither the contract nor candidate arithmetic changed after the first measurement.', '',
        'Next iterations must rerun this entire matched suite, preserve failures and add fresh '
        'holdouts. Promote a candidate only after calibrated full-pipeline costs improve; '
        'larger-group scaling and runtime confidence intervals remain separate outstanding gates.', '']
    return '\n'.join(text)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with args.output.open('x') as out:
        out.write(render(json.loads(args.comparison.read_text())))
