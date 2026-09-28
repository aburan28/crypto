"""Render the note and static scoreboard section from the frozen result export."""
import argparse
import html
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SCOREBOARD = REPO / 'docs/index-calculus-scoreboard.html'
URL = 'https://github.com/aburan28/crypto/blob/main/research/ic_candidate_tournament_20260915/goal_20260924/generic-reference-qualification/'
BEGIN = '<!-- BEGIN ic-generic-reference-qualification-20260926 -->'
END = '<!-- END ic-generic-reference-qualification-20260926 -->'


def display(value):
    if value is None:
        return 'unknown / not applicable'
    if isinstance(value, list):
        return '[' + ', '.join(display(v) for v in value) + ']'
    return format(value, '.7g') if isinstance(value, float) else str(value)


def sections(data):
    rows = data['summary']
    result = []
    def prose(title, text):
        result.append((title, text, None))
    def table(title, columns, values):
        result.append((title, '', (columns, values)))
    prose('Measured development comparison; no promoted winner',
        'The registered run 36290704597 completed all 1,350 native/profile pairs '
        '(30 A/A, 330 smoke and 990 development), with zero research failures. '
        'The whole-mode observer study completed all 360 pairs, with zero incomplete or failed pairs. '
        'This is accounting and reference qualification on public synthetic toy groups. '
        'It consumes no improvement round and does not change the accepted reference binding. '
        'One of three improvement rounds is complete; two remain. No global optimum or held-out improvement is established.')
    prose('Decision and comparison scope',
        'The optimized pairinv incumbent remains the cold IC leader. Prepared both is the online '
        'point-estimate leader, at 0.9770794 times incumbent time, with a descriptive 95% interval '
        '[0.9176026, 1.043898]; this does not establish an IC improvement. Generic dense and sparse '
        'IC are roughly 3.7 times slower online, ten times slower in cold native time and 37 times '
        'more costly in cold instructions than the optimized incumbent. Future rounds must retain '
        'the optimized reference. Choosing the generic worker as a weaker baseline would not establish progress. '
        'The frozen selectors choose rho_incumbent_8 for cold instructions and rho_generic_dense_16 '
        'for online time. These are development selections; a reviewed, versioned binding is still required.')
    prose('Primary metric and uncertainty',
        'Every job solves one supplied public point. Native online time begins after reusable '
        'factor-base, index and log preparation, includes all target-dependent attempts, and ends '
        'after scalar replay. Fixture construction, launch and input loading remain outside this '
        'interval. Three process repetitions are medianed per point; equal-cell geometric means '
        'combine 15 development points across five cells. The repetitions are not 45 independent '
        'targets and this is not shared-target amortization. Intervals are descriptive paired '
        'development bootstraps, with the frozen estimator and seed; they do not correct for '
        'selecting leaders on this panel. IC/rho intervals resample their paired observations '
        'directly, rather than dividing marginal confidence endpoints.')
    table('Single-target online wall time: all 22 configurations',
        ['Variant', 'Online ms', '/ incumbent', 'Descriptive 95% interval',
         'Selected online rho / variant', 'Paired rho/IC 95% interval', 'Verified / scheduled'],
        [[r['alias'], r['online_ms'], r['online_over_incumbent'], r['online_ratio_ci95'],
          r['selected_rho_over_online'], r['selected_rho_over_online_ci95'],
          f"{r['verified_runs']}/{r['scheduled_runs']}"] for r in rows])
    prose('Paired supplied points',
        'The next table shows the selected online IC and rho on each of the 15 development points. '
        'Every interval is one supplied point after reusable preparation through scalar replay. '
        'Full IC1 candidate IDs, RHO1 reference IDs, workload IDs, execution IDs and all four IC '
        'comparisons are retained in RESULTS.json; aliases below map to those exact manifests. '
        'Speedup is rho_online_ms / IC_online_ms, limited to this charged instrumented-worker comparison.')
    table('Selected online IC and rho: point-level observations',
        ['Case', 'Public Q (x,y)', 'Workload', 'IC alias', 'IC ms', 'Rho alias', 'Rho ms', 'Rho / IC', 'Verified'],
        [[r['case'], ','.join(r['public_target']), ','.join(r['workload_ids']), r['arm'],
          r['IC_online_ms'], r['rho_alias'], r['rho_online_ms'], r['online_speedup'], r['verified']]
         for r in data['single_target_online'] if r['arm'] == data['selected']['selected_ic_online']])
    prose('Supplementary complete cold accounting',
        'Cold native time is the full cold child process, including launch, input and reporting tail. '
        'It is not the primary single-target online metric. Cold Ir is Valgrind 3.22 amd64 user-space '
        'guest instructions with exclusive complete phase closure; it is not curve additions. '
        'S = Ir / sqrt(r) uses subgroup order r, not field size. The existing K-instruction floor '
        'applies only to these full-rank relation collectors and is not a universal IC bound; it is '
        'inapplicable to rho. No observer or witness overhead is subtracted.')
    table('Supplementary cold native wall time',
        ['Variant', 'Cold process ms', '/ incumbent', 'Descriptive 95% interval', 'Verified / scheduled'],
        [[r['alias'], r['cold_ms'], r['cold_time_over_incumbent'], r['cold_time_ratio_ci95'],
          f"{r['verified_runs']}/{r['scheduled_runs']}"] for r in rows])
    table('Supplementary complete cold instructions and boundaries',
        ['Variant', 'Cold Ir', 'S = Ir/sqrt(r)', '/ incumbent', 'Descriptive 95% interval',
         '/ selected cold rho', '/ K-instruction floor', 'Verified / scheduled'],
        [[r['alias'], r['cold_Ir'], r['S_Ir_per_sqrt_r'], r['Ir_over_incumbent'], r['Ir_ratio_ci95'],
          r['Ir_over_selected_rho'], r['Ir_over_K_floor'],
          f"{r['verified_runs']}/{r['scheduled_runs']}"] for r in rows])
    prose('Actual factor bases and clipped rho widths',
        'All four IC arms have the following independently checked usable base counts B before '
        'sign/Frobenius folding and effective column counts K. The requested points = 6n is a '
        'construction parameter, never B. Factor-base policy was explicitly a comparison variable. '
        'The raw receipts retain rank, verified yield, unsuccessful PDP attempts, query histories, '
        'matrix work, descent and certificates. Larger unresolved PDP queries remain unresolved; '
        'they do not become unsatisfiability proofs.')
    table('Actual base inventory, shared counts across four IC arms',
        ['Cell', 'Usable points B', 'Folded columns K'],
        [[cell, v['usable_points'], v['folded_columns']] for cell, v in rows[0]['base_inventory'].items()])
    table('Effective rho widths, in n17a1 / n19a0 / n23a0 / n23a1 / n31a0 order',
        ['Variant', 'Observed widths by cell'],
        [[r['alias'], '; '.join(f"{cell}: {display(width)}" for cell, width in r['effective_rho_widths'].items())]
         for r in data['qualification']['table'] if r['mode'] == 'rho'])
    prose('Interpretation of width comparisons',
        'Requested widths 8, 16 and 32 clip to the same widths 1, 1, 4, 5 and 2 on this panel. '
        'They are retained as scheduled configurations, but do not constitute independent algorithms '
        'or additional target samples. A selected width cannot support a general batching advantage.')
    ic = [r for r in rows if r['mode'] == 'ic']
    table('Cold instruction phase shares: stage diagnostics only',
        ['Phase'] + [r['alias'] + ' (%)' for r in ic],
        [[phase] + [format(100 * float(r['cold_instruction_phase_shares'][phase]), '.5f') for r in ic]
         for phase in ic[0]['cold_instruction_phase_shares']])
    prose('What the phase ledger supports',
        'Shares are means of per-run fractions of complete cold Ir, with equal target and cell weight. '
        'Generic precomputation accounts for about 77%; final relation LA accounts for less than '
        '0.2%. The optimized incumbent spends about 1.6% in final relation LA. Eliminating that '
        'stage alone therefore cannot supply a 20% complete-cold improvement on these measurements. '
        'This is prioritization evidence, not attribution of cost to timers or witness generation. '
        'F4/F5 internal matrix reduction and final relation LA remain distinct stages.')
    prose('Whole-mode observer effects',
        'The registered enabled/legacy study checks semantic agreement and retains both outputs. '
        'It changes more than timers: collection mode, tracing and matrix/batch reporting also differ. '
        'Both modes retain witness reporting. Legacy has OBS1 observation identifiers and null '
        'scientific admission and phase costs. The common outer interval omits enabled phase-boundary '
        'bookkeeping, so it cannot replace the fully charged scientific online interval. '
        'Ratios near one do not establish zero or low timer overhead. No overhead is subtracted. '
        'The following 95% intervals are descriptive, from 10,000 paired within-cell resamples; '
        'they do not supply promotion evidence.')
    for metric, name in (('outer_online', 'Common supplied-point outer wall interval'),
                         ('process_wall', 'Whole-process wall time')):
        table('Observer: ' + name,
            ['Variant', 'Enabled / legacy', 'Descriptive 95% interval', 'Complete pairs'],
            [[alias, effect['metrics'][metric]['enabled_over_legacy'],
              effect['metrics'][metric]['descriptive_percentile_95'], 45]
             for alias, effect in data['observer']['effects'].items()])
    aa = data['aa']['comparisons'][0]
    prose('Noise control and environment',
        'The 30-run A/A instruction gate passed. A/A online ratio is '
        + display(aa['online']['candidate_over_baseline']) + ', interval '
        + display(aa['online']['ci95']) + '; cold process ratio is '
        + display(aa['native_wall_candidate_over_baseline']) + ', interval '
        + display(aa['native_wall_ci95']) + '. This shows material wall-time noise. '
        'The experiment used the same Linux x86-64 Azure CI host, Rust 1.94.1, Valgrind 3.22.0, '
        'CPU affinity 3, one Rayon thread, 8 GiB and a 180-second child cap. The generic worker '
        'reported pclmulqdq dispatch. This is a virtualized CI environment; no bare-metal, ARM, '
        'GPU or FPGA result is implied. Host/runtime/build records are retained. The host manifest '
        'does not record the physical CPU model or concurrent host load; neither is inferred.')
    prose('Retained failures, portability and delivery gate',
        'Research failures are zero, but the prerequisite controls retain intentional incomplete '
        'runs, including 12 incomplete mixed-adapter profiles and nine incomplete observer-control '
        'pairs. Nine prerequisite audits and the 360-pair observer replay passed locally without '
        'workers. Strict macOS replay of the main comparison failed on a one-ULP statistic: '
        'development/rho_comparisons/rho_prepared_both_8/per_cell/n23a0 was stored as '
        '1.1843664724946852 and recomputed as 1.1843664724946854. Original Linux bytes and this '
        'failure are preserved. The dedicated Linux evidence CI must freshly restore the archive, '
        'reconstruct every export and pass all eleven frozen audits without worker execution. '
        'The verifier has not been weakened to accommodate the local difference.')
    prose('Next bounded research step',
        'Do not redispatch this qualification or sealed round one. Bind the reviewed observer '
        'evidence and separate cold/online IC and rho leaders in a new versioned reference contract. '
        'Then freeze round two with seed 2026092552 and fresh candidate combinations, retaining '
        'diverse alternatives and an exploration slot rather than discarding every local loser. '
        'Exclude all 25 exposed points in fixtures.json, including A/A and smoke, in addition to '
        'all earlier exposed panels. The n29a1 holdout was not generated or inspected here. '
        'Promotion still requires fresh confirmation and replay under the predeclared familywise '
        'rule, at least 20% lower complete cold Ir and cold native time, no online regression, '
        'no cell regression above 10%, and independently verified answers. Two rounds remain.')
    return result


def render_markdown(parts):
    out = ['# Generic/reference comparison and observer evidence\n',
           '[Frozen protocol](PROTOCOL.md) · [Full result export](RESULTS.json) · '
           '[Exposed fixtures](fixtures.json) · '
           '[Completed workflow](https://github.com/aburan28/crypto/actions/runs/36290704597)\n']
    for title, prose, table in parts:
        out.append('## ' + title + '\n')
        if prose:
            out.append(prose + '\n')
        if table:
            columns, rows = table
            out.append('| ' + ' | '.join(columns) + ' |')
            out.append('| ' + ' | '.join('---' for _ in columns) + ' |')
            out.extend('| ' + ' | '.join(display(v) for v in row) + ' |' for row in rows)
            out.append('')
    entry = next(e for e in json.loads((HERE.parents[1] / 'evidence/manifest.json').read_text())['archives']
                 if e['file'] == 'ic-generic-reference-qualification-20260926.tar.zst')
    out.extend(['## Durable archive and reproduction\n',
        'The archive is committed in this repository, with its identity in '
        '[the evidence manifest](../../evidence/manifest.json). It contains the comparison, '
        'observer study, frozen evaluators, raw profiles, source and controlled build, prerequisite '
        'artifacts, GitHub artifact provenance, workflow record and local audit receipts. '
        'Compression preserves original profile bytes; no measurement is rewritten.\n',
        f"Archive: {entry['file']}; {entry['bytes']:,} bytes, {entry['files']:,} retained files, "
        f"{entry['uncompressed_file_bytes']:,} reconstructed file bytes.\n",
        f"SHA-256: `{entry['sha256']}`.\n",
        'From the repository root, using Python 3.12 and zstd:\n',
        '```sh\npython3.12 -m unittest discover \\\n'
        '  -s research/ic_candidate_tournament_20260915/goal_20260924/generic-reference-qualification \\\n'
        '  -p test_archive.py -v\n```\n',
        'Run on Linux for the complete strict replay; the Linux-only audit test is explicitly '
        'skipped on macOS. No test executes a research worker. The separate evidence workflow '
        'avoids repeating the large historical replay in every producer job.\n',
        'The immutable fixture export has SHA-256 '
        '`b881798bfc56acdd8b8cc52a14b7501572c181a4e25e6bab1c8e7e61a742c7eb`. '
        'RESULTS.json retains all 990 canonical development run records, the original qualification '
        'and observer summaries, exact identities and 60 paired IC/rho point rows. The raw archive '
        'retains the complete 1,350-slot comparison and 360 observer pairs.\n',
        'Regenerate this note and the canonical scoreboard with `render_report.py`; '
        '`render_report.py --check` verifies that both remain exact renderings of the exported evidence.\n'])
    return '\n'.join(out)


def render_html(parts):
    out = [BEGIN, '<section class="panel" id="ic-generic-reference-qualification-20260926">',
           '<h2>Generic/reference qualification: complete evidence, no promoted IC winner</h2>',
           '<p><span class="chip">accounting</span> '
           f'<a href="{URL}README.md">Research note and archive</a> · '
           f'<a href="{URL}RESULTS.json">Frozen results</a> · '
           f'<a href="{URL}PROTOCOL.md">Registered protocol</a></p>']
    for title, prose, table in parts:
        out.append('<h3>' + html.escape(title) + '</h3>')
        if prose:
            out.append('<p>' + html.escape(prose) + '</p>')
        if table:
            columns, rows = table
            out.append('<div class="table-scroll"><table><thead><tr>' +
                       ''.join('<th>' + html.escape(c) + '</th>' for c in columns) +
                       '</tr></thead><tbody>')
            out.extend('<tr>' + ''.join('<td>' + html.escape(display(v)) + '</td>' for v in row) + '</tr>'
                       for row in rows)
            out.append('</tbody></table></div>')
    return '\n'.join(out + ['</section>', END])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    parts = sections(json.loads((HERE / 'RESULTS.json').read_text()))
    note, section = render_markdown(parts), render_html(parts)
    original = SCOREBOARD.read_text()
    if BEGIN in original:
        start = original.index(BEGIN)
        end = original.index(END, start) + len(END)
        page = original[:start] + section + original[end:]
    else:
        marker = '<section class="panel" id="ic-generic-adapter-control-20260926">'
        if original.count(marker) != 1:
            raise ValueError('missing or ambiguous scoreboard insertion point')
        page = original.replace(marker, section + '\n\n' + marker, 1)
    if args.check:
        if (HERE / 'README.md').read_text() != note or original != page:
            raise ValueError('note or scoreboard differs from frozen report')
    else:
        (HERE / 'README.md').write_text(note)
        SCOREBOARD.write_text(page)
    print('VERIFIED' if args.check else 'RENDERED')


if __name__ == '__main__':
    main()
