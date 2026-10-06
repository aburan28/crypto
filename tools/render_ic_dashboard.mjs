#!/usr/bin/env node
// Render existing evidence. No arithmetic research, solver execution or statistics.
import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import { resolve, dirname } from 'node:path';
import assert from 'node:assert/strict';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const PAGE = resolve(ROOT, 'docs/index-calculus-scoreboard.html');
const GOAL = 'research/ic_candidate_tournament_20260915/goal_20260924/';
const START = '<!-- BEGIN IC DASHBOARD OVERVIEW -->';
const END = '<!-- END IC DASHBOARD OVERVIEW -->';
const CSS_START = '/* BEGIN IC OVERVIEW STYLES */';
const CSS_END = '/* END IC OVERVIEW STYLES */';
const sha = bytes => createHash('sha256').update(bytes).digest('hex');
const escape = value => String(value).replace(/[&<>"']/g, ch => ({ '&':'&amp;', '<':'&lt;', '>':'&gt;', '"':'&quot;', "'":'&#39;' }[ch]));
const link = (path, label) => `<a href="https://github.com/aburan28/crypto/blob/main/${GOAL}${path}">${label}</a>`;
function source(path, expected) {
  const bytes = readFileSync(resolve(ROOT, GOAL, path));
  assert.equal(sha(bytes), expected, `Evidence changed: ${path}. Review its claim before refreshing its pin.`);
  return [JSON.parse(bytes), { path: GOAL + path, sha256: expected }];
}
function between(text, start, end) {
  assert.equal(text.split(start).length, 2, `Expected one ${start}`);
  const i = text.indexOf(start) + start.length;
  assert.ok(text.indexOf(end, i) >= i, `Missing ${end}`);
  return text.slice(i, text.indexOf(end, i));
}
function replace(text, start, end, value) {
  between(text, start, end);
  const i = text.indexOf(start) + start.length;
  return text.slice(0, i) + value + text.slice(text.indexOf(end, i));
}
function section(text, id) {
  const opening = [...text.matchAll(/<section\b[^>]*>/g)].find(match => match[0].includes(`id="${id}"`));
  assert.ok(opening, `Missing preserved section ${id}`);
  let depth = 0;
  for (const match of text.slice(opening.index).matchAll(/<\/?section\b[^>]*>/g)) {
    depth += match[0].startsWith('</') ? -1 : 1;
    if (!depth) return text.slice(opening.index, opening.index + match.index + match[0].length);
  }
  assert.fail(`Unclosed section ${id}`);
}

const [round3, roundPin] = source('improvement/round3/RESULTS.json', '8e6b506517c945a6180da9c7e7eb85bef4f6a145948deec06bc4a31093aff143');
const [f5, f5Pin] = source('prepared-f5-v3-control-v1/TERMINAL.json', 'a06b68139abad2934fab8127d2ee69ea126ab09a7679ba6e775c9352db09ba1e');
const [nativeF5, nativeF5Pin] = source('native-f5-control-registration-v1/result-v1/data/independent-audit.json', '994ab5bf1f3153eb6d22425eb15647a152cc1ebd14b68347147c279479a9ed69');
const [oldSat, oldSatPin] = source('prepared-one-target-controls-v1/outcome/result.json', '9def5e0e4a3c3fd2e47252d08e3509a88205691c39a77912cd315428b7bed2eb');
const [sat, satPin] = source('native-sat-control-registration-v1/result-v1/data/independent-audit.json', '3570722341b8ca2cd8794d51829cac1a1f3b8c1ff02a07a6af9377a36d2cbab7');
const [ordinary, ordinaryPin] = source('native-ordinary-registration-v1/result-v1/stage-summary.json', '4a2689bedfdb5443c72239067b7a2997f30c20794c1e774db6b76019892e4d84');
const [satMillion, satMillionPin] = source('native-sat-million-registration-v1/result-v1/original-audit.json', 'e2e9caa59ffa61713e4ffd47d98a05850fe4e22bef9594d9431c2c7b24d283be');
const [f5Target, f5TargetPin] = source('native-f5-one-target-registration-v1/result-v1/audit-original.json', '6754ae85061d60d240cd3d8cb50efd8649354255533fc3b7f49ce654486aeb31');
const [pilotNegative, pilotNegativePin] = source('native-sat-budget-pilot-v1/result-v1/query-000/result.json', '51ca35d7d0afe9a3bf5e41337f325fc29cc81c94f8f3f1205962a1b96635f387');
const [pilotFeasible, pilotFeasiblePin] = source('native-sat-budget-pilot-v1/result-v1/query-015/result.json', 'b4bd23aa4a5d6d0f277a80c5d5f011041918d6812adc390081debc0b6865552d');
assert.equal(ordinary.stage_only, true);
assert.equal(ordinary.ordinary_queries, 512);
assert.equal(ordinary.f5.final_rank, 29);
assert.equal(ordinary.f5.folded_columns, 29);
assert.equal(ordinary.f5.independently_recovered_logs, 29);
assert.deepEqual(ordinary.f5.outcome_mix, {proved_unsat:383,witness:129});
assert.equal(ordinary.cms.final_rank, 22);
assert.equal(ordinary.cms.folded_columns, 29);
assert.equal(ordinary.cms.independently_recovered_logs, null);
assert.deepEqual(ordinary.cms.outcome_mix, {incomplete:490,witness:22});
assert.equal(ordinary.cms.feasible_inconclusive_queries, 107);
assert.equal(ordinary.online_speedup, null);
assert.equal(satMillion.status, 'PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT');
assert.equal(satMillion.source_bound_execution_admitted, true);
assert.equal(satMillion.audited_ordinary_queries, 512);
assert.equal(satMillion.mathematics.mathematical_preparation_complete, true);
assert.equal(satMillion.mathematics.rank, 29);
assert.equal(satMillion.mathematics.folded_columns, 29);
assert.equal(satMillion.mathematics.usable_points, 62);
assert.equal(satMillion.mathematics.independently_recovered_column_logs.length, 29);
assert.deepEqual(satMillion.mathematics.ordinary_outcome_mix,
  {incomplete:15,proved_unsat:339,timeout:13,witness:145});
assert.equal(satMillion.mathematics.accepted_rows, 145);
assert.equal(satMillion.online_speedup, null);
assert.equal(f5Target.status, 'PASS_NATIVE_SOURCE_BOUND_TARGET_EXECUTION_AUDIT');
assert.equal(f5Target.source_bound_execution_admitted, true);
assert.equal(f5Target.verified_recovery, true);
assert.equal(f5Target.recorded_online_interval_ns, 6501300958);
assert.equal(f5Target.headline_eligible, false);
assert.equal(f5Target.online_speedup, null);
assert.equal(pilotNegative.cold.status, 'SOURCE_UNSAT');
assert.equal(pilotNegative.stdin.status, 'SOURCE_UNSAT');
assert.equal(pilotFeasible.cold.status, 'SAT_MODEL');
assert.equal(pilotFeasible.stdin.status, 'SAT_MODEL');
assert.equal(pilotNegative.prestarted_stdin_compatible, true);
assert.equal(pilotFeasible.prestarted_stdin_compatible, true);
assert.equal(pilotNegative.source_bound_target_admitted, false);
assert.equal(pilotFeasible.source_bound_target_admitted, false);
const solverComparisonPath = 'research/ic_solver_online_20261003/confirmation/RESULT.md';
const solverComparisonSha256 = 'e4cfed0199e4d3c55731dfcd659a160d0af73773f13f19d06e0813bdfd123ffe';
assert.equal(sha(readFileSync(resolve(ROOT, solverComparisonPath))), solverComparisonSha256,
  'Fresh solver comparison changed; review its result before refreshing the dashboard');
const solverComparisonPin = { path: solverComparisonPath, sha256: solverComparisonSha256 };
const sharedRankPath = 'research/notes/ecc2k130/shared_rank_folded_table_20261003/EVIDENCE.json';
const sharedRankBytes = readFileSync(resolve(ROOT, sharedRankPath));
const sharedRankSha = 'cf61ec46ae7be8eb58559a17a6f90eb91cfc1f1b67dae72beaabcebed6ba2531';
assert.equal(sha(sharedRankBytes), sharedRankSha, 'Shared-rank evidence changed; review the claim before refreshing its pin.');
const sharedRank = JSON.parse(sharedRankBytes);
assert.equal(sharedRank.decision, 'TARGET_BLIND_RANK42_ADMITTED_N37_SOURCE');
assert.equal(sharedRank.class, 'stage_diagnostic');
assert.equal(sharedRank.counts.trials, 55);
assert.equal(sharedRank.counts.rank, 42);
assert.equal(sharedRank.counts.column_logs_verified, 42);
assert.deepEqual(sharedRank.gae_lower_bound, {
  base: 4260, table: 178710, rank_search: 8498,
  linear_algebra: 1110, verification: 2248, total: 194826
});
assert.equal(sharedRank.rho_ratio, null);
const sharedRankPin = {path: sharedRankPath, sha256: sharedRankSha};
const sharedTargetsPath = 'research/notes/ecc2k130/n37_shared_rank_targets_20261003/EVIDENCE.json';
const sharedTargetsBytes = readFileSync(resolve(ROOT, sharedTargetsPath));
const sharedTargetsSha = '7b8f9c0d3fce0da4165d482aeaea3d4698c634173df53efaa936f17ce8f25c27';
assert.equal(sha(sharedTargetsBytes), sharedTargetsSha, 'Shared-target evidence changed; review the claim before refreshing its pin.');
const sharedTargets = JSON.parse(sharedTargetsBytes);
assert.equal(sharedTargets.decision, 'SHARED_RANK_TARGET16_RECOVERY_ADMITTED_N37_SOURCE');
assert.equal(sharedTargets.class, 'accounting');
assert.equal(sharedTargets.counts.targets, 16);
assert.equal(sharedTargets.counts.verified_scalar_recoveries, 16);
assert.equal(sharedTargets.counts.direct_three_summand_hits, 16);
assert.equal(sharedTargets.gae_lower_bound.cold_total, 197454);
assert.equal(sharedTargets.rho_ratio, null);
const sharedTargetsPin = {path: sharedTargetsPath, sha256: sharedTargetsSha};
const sharedEcbenchPath = 'research/ecbench_n37_shared_rank_20261004/RESULT_FINAL.json';
const sharedEcbenchBytes = readFileSync(resolve(ROOT, sharedEcbenchPath));
const sharedEcbenchSha = 'cabeb087e0991580831434ba35c688ace9072e346cd7ceeb0522aeb5c0578ede';
assert.equal(sha(sharedEcbenchBytes), sharedEcbenchSha, 'One-target ecbench result changed; review the decision before refreshing its pin.');
const sharedEcbench = JSON.parse(sharedEcbenchBytes);
assert.equal(sharedEcbench.status, 'verified_bounded_l0_diagnostic');
assert.equal(sharedEcbench.all_measured_verified, true);
assert.equal(sharedEcbench.target_count, 1);
assert.equal(sharedEcbench.measured_rounds, 5);
assert.equal(sharedEcbench.online_speedup, null);
assert.equal(sharedEcbench.cold_counted_ic_over_rho_bounded, 13.590881825669284);
const sharedEcbenchPin = {path: sharedEcbenchPath, sha256: sharedEcbenchSha};
const sharedComparisonPath = 'research/ecbench_n37_shared_rank_20261004/sessions/mac_arm64_l0_02/comparisons/rho-strong__ic-shared.json';
const sharedComparisonBytes = readFileSync(resolve(ROOT, sharedComparisonPath));
const sharedComparisonSha = '28f375249d4277b44ebc1920691acd0b242dbf8a1c4dd8a69ec352b31d4fe810';
assert.equal(sha(sharedComparisonBytes), sharedComparisonSha, 'Shared-rank comparison changed; review before refreshing its pin.');
const sharedComparison = JSON.parse(sharedComparisonBytes);
assert.equal(sharedComparison.ops.bounded, true);
assert.equal(sharedComparison.ops.pairs, 5);
assert.equal(sharedComparison.wall.status, 'descriptive');
assert.equal(sharedComparison.ops.ratio_b_over_a, sharedEcbench.cold_counted_ic_over_rho_bounded);
const sharedComparisonPin = {path: sharedComparisonPath, sha256: sharedComparisonSha};
const sharedControlPath = 'research/ecbench_n37_shared_rank_20261004/sessions/mac_arm64_l0_02/comparisons/ic-shared__ic-control.json';
const sharedControlBytes = readFileSync(resolve(ROOT, sharedControlPath));
const sharedControlSha = 'd5277601068a516d72131ebef5f024d2d305cfbe577cd577ce7b06370742573d';
assert.equal(sha(sharedControlBytes), sharedControlSha, 'Shared-rank A/A comparison changed; review before refreshing its pin.');
const sharedControl = JSON.parse(sharedControlBytes);
assert.equal(sharedControl.ops.ratio_b_over_a, 1);
const sharedControlPin = {path: sharedControlPath, sha256: sharedControlSha};
const sharedIndependentPath = 'research/ecbench_n37_shared_rank_20261004/independent_validation_20261004/RECEIPT.json';
const sharedIndependentBytes = readFileSync(resolve(ROOT, sharedIndependentPath));
const sharedIndependentSha = 'bc2e2ebf42f248d55d81691ffb87d337eb0063c3c8e9bd516bd654bbde56729f';
assert.equal(sha(sharedIndependentBytes), sharedIndependentSha, 'Shared-rank independent replay changed; review its claim before refreshing the dashboard.');
const sharedIndependent = JSON.parse(sharedIndependentBytes);
assert.equal(sharedIndependent.ok, true);
assert.equal(sharedIndependent.session_id, 'ECBS1h33c3e5e4f75d');
assert.equal(sharedIndependent.auditor_env_class_id, 'ECBENV2hda64b23436f4');
assert.equal(sharedIndependent.replays.length, 15);
assert.ok(sharedIndependent.replays.every(row => row.reproduced));
assert.equal(sharedIndependent.files['records.jsonl'], '8046061f5078f254627a343998b810c7e570fc9fe575f403466114e70a0b1693');
const sharedIndependentPin = {path: sharedIndependentPath, sha256: sharedIndependentSha};
const sharedIndependentClaimPath = 'research/ecbench_n37_shared_rank_20261004/independent_validation_20261004/CLAIM_INDEPENDENT_DIAGNOSTIC.json';
const sharedIndependentClaimBytes = readFileSync(resolve(ROOT, sharedIndependentClaimPath));
const sharedIndependentClaimSha = 'adabcd19c30e81c6f3967377a13e8ef2acd23ade02eb26ca77b5fa2312d8ab5c';
assert.equal(sha(sharedIndependentClaimBytes), sharedIndependentClaimSha, 'Shared-rank attached claim changed; review before refreshing the dashboard.');
const sharedIndependentClaim = JSON.parse(sharedIndependentClaimBytes);
assert.equal(sharedIndependentClaim.independent_validation, true);
assert.equal(sharedIndependentClaim.independent_replay.receipt_sha256, sharedIndependentSha);
assert.equal(sharedIndependentClaim.isolation_levels.ic, 'L0');
assert.equal(sharedIndependentClaim.isolation_levels.rho, 'L0');
assert.equal(sharedEcbench.online_speedup, null);
const sharedIndependentClaimPin = {path:sharedIndependentClaimPath, sha256:sharedIndependentClaimSha};
const sharedIndependentCheckPath = 'research/ecbench_n37_shared_rank_20261004/independent_validation_20261004/CHECK.json';
const sharedIndependentCheckBytes = readFileSync(resolve(ROOT, sharedIndependentCheckPath));
const sharedIndependentCheckSha = '06a4356519163ec6d8d62f90fdce64fcba509693fd924cc98ab00dec80e3e482';
assert.equal(sha(sharedIndependentCheckBytes), sharedIndependentCheckSha, 'Shared-rank independent claim check changed; review before refreshing the dashboard.');
const sharedIndependentCheck = JSON.parse(sharedIndependentCheckBytes);
assert.equal(sharedIndependentCheck.status, 'PASS');
const sharedIndependentCheckPin = {path:sharedIndependentCheckPath, sha256:sharedIndependentCheckSha};
const rankColumnsPath = 'research/ecbench_n37_rank_columns_20261004/RESULT.json';
const rankColumnsBytes = readFileSync(resolve(ROOT, rankColumnsPath));
const rankColumnsSha = '9d3ba4708f70e2796bbdb8a541a85b794972563d6ecfe70bde160a698da49fbb';
assert.equal(sha(rankColumnsBytes), rankColumnsSha, 'n37 column decision changed; review its evidence before refreshing the dashboard.');
const rankColumns = JSON.parse(rankColumnsBytes);
assert.equal(rankColumns.status, 'independently_replayed_l0_bounded_diagnostic');
assert.equal(rankColumns.selection.decision, 'COUNTED_ENGINEERING_LEAD');
assert.equal(rankColumns.selection.selected_k, 16);
assert.equal(rankColumns.target_count, 8);
assert.equal(rankColumns.measured_rounds_per_target, 5);
assert.equal(rankColumns.candidates.length, 6);
assert.equal(rankColumns.online_speedup, null);
assert.equal(rankColumns.fully_priced_cold_speedup, null);
const rankColumnsPin = {path:rankColumnsPath, sha256:rankColumnsSha};
const rankColumnsAuditPath = 'research/ecbench_n37_rank_columns_20261004/independent_validation/RECEIPT.json';
const rankColumnsAuditBytes = readFileSync(resolve(ROOT, rankColumnsAuditPath));
const rankColumnsAuditSha = 'efb14a5cb1c60291accc4d28518178f6a7d7ed5305178459936252b346b94029';
assert.equal(sha(rankColumnsAuditBytes), rankColumnsAuditSha, 'n37 column independent receipt changed.');
const rankColumnsAudit = JSON.parse(rankColumnsAuditBytes);
assert.equal(rankColumnsAudit.ok, true);
assert.equal(rankColumnsAudit.session_id, rankColumns.session_id);
assert.equal(rankColumnsAudit.replays.length, 320);
assert.ok(rankColumnsAudit.replays.every(row => row.reproduced));
assert.equal(rankColumns.independent_audit_sha256, rankColumnsAuditSha);
const rankColumnsAuditPin = {path:rankColumnsAuditPath, sha256:rankColumnsAuditSha};
const callgrindPath = 'research/ecbench_callgrind_solve_20261004/DECISION.json';
const callgrindBytes = readFileSync(resolve(ROOT, callgrindPath));
const callgrindSha = 'c4bd4b3b8e1e95520bca5507706b2434eb8bdfb7d7281fdfd6998a71b18dd3aa';
assert.equal(sha(callgrindBytes), callgrindSha, 'n37 whole-solve instruction decision changed.');
const callgrind = JSON.parse(callgrindBytes);
assert.equal(callgrind.unit, 'callgrind.Ir');
assert.equal(callgrind.workloads, 8);
assert.equal(callgrind.profiles, 32);
assert.equal(callgrind.all_archived_and_profiled_scalars_verified, true);
assert.equal(callgrind.decision, 'retain_k16_instruction_lead');
assert.equal(callgrind.comparisons[0].numerator, 'ic-k16');
assert.equal(callgrind.comparisons[0].denominator, 'ic-k42');
assert.equal(callgrind.comparisons[1].denominator, 'rho-strong');
const callgrindPin = {path:callgrindPath, sha256:callgrindSha};
const k8ConfirmPath = 'research/ecbench_n37_k8_k16_20261004/DECISION.json';
const k8ConfirmBytes = readFileSync(resolve(ROOT, k8ConfirmPath));
const k8ConfirmSha = 'c4a33ffd9045a479ffe2f2a0ae2ecb60519cb293ac326d4ac84b4014f7980127';
assert.equal(sha(k8ConfirmBytes), k8ConfirmSha, 'n37 K8/K16 instruction decision changed.');
const k8Confirm = JSON.parse(k8ConfirmBytes);
assert.equal(k8Confirm.schema, 'ecbench.k8_k16_callgrind_decision/v1');
assert.equal(k8Confirm.unit, 'callgrind.Ir');
assert.equal(k8Confirm.workloads, 16);
assert.equal(k8Confirm.profiles, 64);
assert.equal(k8Confirm.all_archived_and_profiled_scalars_verified, true);
assert.equal(k8Confirm.decision, 'select_k8_for_larger_field_gate');
assert.equal(k8Confirm.comparisons[0].numerator, 'ic-k16');
assert.equal(k8Confirm.comparisons[0].denominator, 'ic-k8');
assert.equal(k8Confirm.comparisons[1].denominator, 'rho-strong');
const k8ConfirmPin = {path:k8ConfirmPath, sha256:k8ConfirmSha};
const onlineIrPath = 'research/ecbench_n37_online_ir_20261004/DECISION.json';
const onlineIrBytes = readFileSync(resolve(ROOT, onlineIrPath));
const onlineIrSha = 'e711979f64d08ea3417844b980e54e0fe095b4802bb803cf4f726c555be0e36a';
assert.equal(sha(onlineIrBytes), onlineIrSha, 'n37 target-only instruction decision changed.');
const onlineIr = JSON.parse(onlineIrBytes);
assert.equal(onlineIr.schema, 'ecbench.n37_online_ir_decision/v1');
assert.equal(onlineIr.unit, 'callgrind.Ir');
assert.equal(onlineIr.workloads, 16);
assert.equal(onlineIr.profiles, 64);
assert.equal(onlineIr.all_archived_and_profiled_scalars_verified, true);
assert.equal(onlineIr.decision, 'prioritize_k16_for_isolated_n37_online_wall_gate');
assert.equal(onlineIr.online.comparisons[0].numerator, 'ic-k8');
assert.equal(onlineIr.online.comparisons[0].denominator, 'ic-k16');
assert.ok(onlineIr.online.comparisons[0].bootstrap_95[0] > 1.10);
assert.equal(onlineIr.online.aa_max_absolute_relative_deviation, 0);
const onlineIrPin = {path:onlineIrPath, sha256:onlineIrSha};
const onlineIrReceiptPath = 'research/ecbench_n37_online_ir_20261004/independent_validation/RECEIPT.json';
const onlineIrReceiptBytes = readFileSync(resolve(ROOT, onlineIrReceiptPath));
const onlineIrReceiptSha = 'b75740f1a28dd4a5352c15e0e738a82bf9c8d3402e3f6ea797d129c424e05472';
assert.equal(sha(onlineIrReceiptBytes), onlineIrReceiptSha, 'n37 target-only independent receipt changed.');
const onlineIrReceipt = JSON.parse(onlineIrReceiptBytes);
assert.equal(onlineIrReceipt.ok, true);
assert.equal(onlineIrReceipt.records, 384);
assert.equal(onlineIrReceipt.verified_records, 384);
assert.equal(onlineIrReceipt.replays.length, 320);
assert.ok(onlineIrReceipt.replays.every(row => row.reproduced));
assert.equal(onlineIrReceiptSha, onlineIr.independent_audit_sha256);
const onlineIrReceiptPin = {path:onlineIrReceiptPath, sha256:onlineIrReceiptSha};
const nativeWallPath = 'research/ecbench_n37_native_online_wall_20261004/DECISION.json';
const nativeWallBytes = readFileSync(resolve(ROOT, nativeWallPath));
const nativeWallSha = 'eebbe2eef9524c7f45a5bfd1c0a2b312d72063a3a8865140a6e25200235c21e2';
assert.equal(sha(nativeWallBytes), nativeWallSha, 'n37 native wall decision changed.');
const nativeWall = JSON.parse(nativeWallBytes);
assert.equal(nativeWall.schema, 'ecbench.native_online_wall_decision/v1');
assert.equal(nativeWall.status, 'complete_exploratory_hosted');
assert.equal(nativeWall.decision, 'carry_both_host_noise_exceeds_gate');
assert.equal(nativeWall.measured_verified_runs, 320);
assert.equal(nativeWall.same_target_paired_rounds, 80);
assert.deepEqual(nativeWall.isolation_levels, {L1: 320});
assert.equal(nativeWall.all_measured_L2, false);
assert.equal(nativeWall.online_speedup, null);
assert.ok(nativeWall.k16_aa_max_relative_deviation > .05);
assert.equal(nativeWall.primary_one_target.target_index, 0);
const nativeWallPin = {path:nativeWallPath, sha256:nativeWallSha};
const nativeWallEvidencePath = 'research/ecbench_n37_native_online_wall_20261004/EVIDENCE.json';
const nativeWallEvidenceBytes = readFileSync(resolve(ROOT, nativeWallEvidencePath));
const nativeWallEvidenceSha = 'bb4dd48030029d21019a21db3cd10e646ef815d727db084013b7bdba3ed8d422';
assert.equal(sha(nativeWallEvidenceBytes), nativeWallEvidenceSha, 'n37 native wall evidence changed.');
const nativeWallEvidence = JSON.parse(nativeWallEvidenceBytes);
assert.equal(nativeWallEvidence.decision_sha256, nativeWallSha);
assert.equal(nativeWallEvidence.records_sha256, nativeWall.records_sha256);
assert.equal(nativeWallEvidence.independent_receipt_sha256, nativeWall.independent_receipt_sha256);
assert.equal(sha(readFileSync(resolve(ROOT, 'research/ecbench_n37_native_online_wall_20261004/SHA256SUMS'))), nativeWallEvidence.raw_sha256s_sha256);
assert.equal(sha(readFileSync(resolve(ROOT, 'research/ecbench_n37_native_online_wall_20261004/SOURCE-ARTIFACT-SHA256SUMS'))), nativeWallEvidence.source_artifact_sha256s_sha256);
assert.equal(nativeWallEvidence.admitted_online_speedup, null);
const nativeWallEvidencePin = {path:nativeWallEvidencePath, sha256:nativeWallEvidenceSha};
assert.equal(round3.decision.winner, 'incumbent');
assert.equal(round3.decision.promotion_eligible, false);
assert.ok(f5.scalar_verified && !f5.headline_online_admissible && !f5.fresh_paired_qualification);
assert.equal(oldSat.arms.sat.status, 'ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL');
assert.ok(sat.source_bound_execution_admitted && sat.target_complete && sat.scalar_verified);
assert.ok(!sat.headline_eligible && !sat.fresh_paired_qualification && !sat.promotion_eligible);
assert.equal(sat.online_speedup, null);
assert.equal(sat.audited_attempts.length, 3);
assert.equal(sat.preparation_admission.new_queries, 0);
assert.ok(nativeF5.source_bound_execution_admitted && nativeF5.target_complete);
assert.ok(!nativeF5.headline_eligible && !nativeF5.fresh_paired_qualification && !nativeF5.promotion_eligible && !nativeF5.full_goal_complete);
assert.equal(nativeF5.online_speedup, null);
assert.equal(nativeF5.ordinary_queries_executed, 0);
assert.equal(nativeF5.mathematics.audited_attempts.length, 3);
assert.equal(nativeF5.mathematics.preparation_admission.new_queries, 0);
assert.equal(nativeF5.mathematics.recovered_scalar, 24886);
const confirmation = round3.decision.confirmation;
const online = confirmation.online;
const ratio = online.candidate_over_baseline;
const [lo, hi] = online.ci95;
const upper = confirmation.familywise.metrics.online_ns.upper;
const rho = round3.decision.winner_over_online_rho.confirmation.online;
const rows = Object.fromEntries(round3.stages.confirmation.table.map(row => [row.alias, row]));
for (const alias of ['incumbent', 'stop7_word', 'rho_online']) {
  assert.ok(rows[alias].complete && rows[alias].verified === rows[alias].scheduled);
}
const f = (value, digits = 3) => value.toFixed(digits);
// Coordinates only: all estimates and intervals above are already recorded.
const percent = value => (value - .85) / .25 * 100;
const forest = `<div class="comparison-plot" data-measurement-chart role="group" aria-label="Challenger to reference ratios; the 1.00 reference line divides faster and slower costs"><div class="plot-row plot-head"><span>Curve cell</span><span class="plot-directions"><span>Faster</span><span>Slower</span></span><span>Ratio</span></div>${Object.entries(online.per_cell).map(([cell,value])=>`<div class="plot-row"><span>${escape(cell)}</span><span class="plot-track"><i class="plot-link" style="left:${Math.min(60,percent(value))}%;width:${Math.abs(60-percent(value))}%"></i><i class="plot-dot ${value>1?'plot-slower':'plot-faster'}" style="left:${percent(value)}%"></i></span><strong>${f(value)}×</strong></div>`).join('')}<div class="plot-row plot-overall"><strong>Overall</strong><span class="plot-track"><i class="plot-interval" style="left:${percent(lo)}%;width:${percent(hi)-percent(lo)}%"></i><i class="plot-dot" style="left:${percent(ratio)}%"></i></span><strong>${f(ratio)}×</strong></div><div class="plot-row plot-head"><span></span><span class="plot-directions"><span>0.85×</span><span>1.10×</span></span><span></span></div><p class="plot-reference">Dashed line = 1.00× reference time</p></div>`;
const rhoGraph = `<div class="rho-plot" data-measurement-chart role="group" aria-label="Recorded single-target online milliseconds, shown on a zero-based linear scale">${['incumbent','stop7_word','rho_online'].map((alias,index)=>`<div class="rho-row"><div><strong>${['Incumbent IC','Challenger IC','Matched rho'][index]}</strong><span>${f(rows[alias].online_ms,6)} ms</span></div><div class="rho-track"><span style="width:${rows[alias].online_ms/rows.rho_online.online_ms*100}%;background:${index===2?'var(--lab-muted)':'var(--lab-blue)'}"></span></div></div>`).join('')}<div class="rho-axis"><span>0</span><span>0.253 ms</span></div></div>`;

let page = readFileSync(PAGE, 'utf8');
const oldFront = between(page, START, END);
const historical = section(oldFront, 'lab-best');
const progressRaw = readFileSync(resolve(ROOT,'docs/ic/progress-timeline.json'));
const retainedPanels = [
  ['cold-compact-orbit-20261003','Historical compact-orbit cold-cost panel'],
  ['full-rank-compact-orbit-20261003','Historical compact-orbit full-rank panel'],
  ['n61-compact-orbit-hosted-20261003','n61 compact-orbit autolab panels · hosted isolated reruns'],
  ['lab-ecbench-all','Cross-method evidence · every measured candidate'],
  ['lab-progress','Research progress over time · separate regimes and references'],
  ['lab-currency','How new results update this dashboard']
].map(([id,label])=>{
  let html = section(oldFront,id);
  if (id === 'lab-progress') {
    const marker = /(<script type="application\/json" id="progress-data">)[\s\S]*?(<\/script>)/;
    assert.ok(marker.test(html), 'Progress panel has no canonical data embed');
    html = html.replace(marker, (_, start, end) => start + progressRaw.toString('utf8').trim() + end);
  }
  return {id,label,html};
});
const historicalDetails = retainedPanels.map(panel=>`<details class="dash-details"><summary>${panel.label}</summary><div class="dash-detail-body">${panel.html}</div></details>`).join('\n');
// The legacy collapse handler must not hide the overview's graphs. Its only
// permitted migration is a selector scope change; all evidence stays verbatim.
const oldSelector = 'document.querySelectorAll(SEL)';
const scopedSelector = "document.getElementById('legacy-evidence').querySelectorAll(SEL)";
const ledgerBefore = between(page, '<div id="legacy-evidence">', '<!-- END IC EVIDENCE LIBRARY -->').replace(oldSelector, scopedSelector);
page = page.replace(oldSelector, scopedSelector);

const stagePercent = count => (100 * count / ordinary.ordinary_queries).toFixed(3);
const naturalStageChart = `<figure class="dash-chart stage-figure"><figcaption><strong>What happened on the same 512 ordinary queries?</strong><span>Each bar includes all 512 attempts · different solver limits · target-free preparation only</span></figcaption>
<div class="stage-row"><div><strong>Matrix F5</strong><span>129 verified witnesses · 383 proved negatives</span></div><div class="stage-track"><i class="stage-witness" style="width:${stagePercent(ordinary.f5.outcome_mix.witness)}%"></i><i class="stage-negative" style="width:${stagePercent(ordinary.f5.outcome_mix.proved_unsat)}%"></i></div><b>Rank ${ordinary.f5.final_rank}/${ordinary.f5.folded_columns}</b></div>
<div class="stage-row"><div><strong>CryptoMiniSat · 100k conflicts</strong><span>22 witnesses · 490 budget-inconclusive</span></div><div class="stage-track"><i class="stage-witness" style="width:${stagePercent(ordinary.cms.outcome_mix.witness)}%"></i><i class="stage-incomplete" style="width:${stagePercent(ordinary.cms.outcome_mix.incomplete)}%"></i></div><b>Rank ${ordinary.cms.final_rank}/${ordinary.cms.folded_columns}</b></div>
<div class="stage-row"><div><strong>CryptoMiniSat · 1m conflicts</strong><span>145 witnesses · 339 proved negatives · 15 incomplete · 13 timeouts</span></div><div class="stage-track"><i class="stage-witness" style="width:${stagePercent(satMillion.mathematics.ordinary_outcome_mix.witness)}%"></i><i class="stage-negative" style="width:${stagePercent(satMillion.mathematics.ordinary_outcome_mix.proved_unsat)}%"></i><i class="stage-incomplete" style="width:${stagePercent(satMillion.mathematics.ordinary_outcome_mix.incomplete + satMillion.mathematics.ordinary_outcome_mix.timeout)}%"></i></div><b>Rank ${satMillion.mathematics.rank}/${satMillion.mathematics.folded_columns}</b></div>
<p class="stage-legend"><span><i class="stage-witness"></i> Verified witness</span><span><i class="stage-negative"></i> Proved no three-sum</span><span><i class="stage-incomplete"></i> Unresolved: budget or timeout</span></p>
<p class="chart-note">The frozen 1m-conflict SAT panel independently recovered all 29 logs. Its 28 unresolved attempts stay unresolved; rank 29/29 is a reusable-preparation result, not a target solve. The 100k-conflict panel's shorter incomplete work was never a speed win. Solver limits differ, so these bars show coverage, not comparative speed. ${link('native-sat-million-registration-v1/result-v1/RESULT.md','New SAT audit →')} ${link('native-ordinary-registration-v1/result-v1/RESULT.md','Earlier paired panel →')}</p></figure>`;

const front = `
<main class="dash" id="ic-overview">
<a class="skip-link" href="#lab-results">Skip to measured results</a>
<nav class="dash-nav" aria-label="Dashboard"><a class="dash-brand" href="#ic-overview"><span class="brand-mark">IC</span> Research lab</a><div><a href="#lab-results">Results</a><a href="#lab-readiness">Solver readiness</a><a href="#lab-pipeline">Pipeline</a><a href="#evidence-search">Evidence</a><a href="../browser/">Lab browser</a></div></nav>
<header class="dash-hero"><div><p class="dash-kicker">Index calculus · bounded autolab · evidence snapshot 5 October 2026</p><h1>Is the next candidate<br>actually better?</h1><p class="dash-lead">Compare complete solutions to the same one target. Keep failed attempts, check the answer, and make uncertainty part of the decision.</p></div><aside class="hero-decision" aria-label="Current decision"><span class="status">Current decision</span><strong>Keep the incumbent</strong><p>No challenger passed the last tournament’s promotion gate. F5 and SAT still need a fresh paired comparison.</p><a href="#lab-next">See what comes next →</a></aside></header>
<div class="dash-metrics" aria-label="At a glance"><article><span class="status">Tournament</span><h2>3 rounds closed</h2><p>The incumbent remains. Historical confirmation targets stay closed.</p></article><article><span class="status good">Disclosed-input correctness</span><h2>F5 solved one registered point</h2><p>Its audited natural panel reached rank 29/29, then its frozen controller recovered and verified one public target. This is not a fresh speed ranking.</p></article><article><span class="status good">SAT preparation</span><h2>All 29 logs recovered</h2><p>The new 512-query, 1m-conflict panel reached rank 29/29. A fresh target solve and same-point comparison are still pending.</p></article></div>

<section class="dash-section" id="lab-results" aria-labelledby="results-title"><div class="dash-section-head"><div><p class="dash-kicker">01 / Tournament result</p><h2 id="results-title">A promising average. No reliable win.</h2></div><span class="status pending">Incumbent retained</span></div>
<div class="dash-result-grid"><figure class="dash-chart"><figcaption id="ratio-title"><strong>Last challenger vs qualified IC reference</strong><span id="ratio-desc">Online time ratio · below 1 is faster · one point per solve</span></figcaption><div class="chart-scroll" tabindex="0" aria-label="Scrollable challenger comparison graph">${forest}</div><p class="chart-note">Circles: faster cells. Diamonds: slower cells. Only the overall mark has a descriptive 95% interval; individual dots are estimates.</p></figure><div class="dash-explanation"><span class="result-number">${f(ratio)}× <span>reference time</span></span><h3>The interval crosses “no improvement.”</h3><p>Recorded 95% interval: <strong>${f(lo)}–${f(hi)}×</strong>. Two of six curve cells were slower.</p><p class="decision-note panel-summary"><strong>Decision:</strong> the stricter familywise upper bound is <strong>${f(upper)}×</strong>. The frozen promotion gate failed.</p><p><strong>${confirmation.paired_cases} paired targets · 6 curve cells</strong><br>Round 3 confirmation, synthetic toy panel. The challenger is <code>stop7_word</code>; its denominator is the qualified <code>ic_online</code> role.</p>${link('improvement/round3/README.md','Read the decision and full table →')}</div></div>
<details class="dash-details"><summary>Exact ratios, timing boundary and source evidence</summary><div class="dash-detail-body"><p>The online interval begins after reusable preparation and ends after scalar replay. Cell names indicate field degree/model, not subgroup bits. These are IC implementation comparisons, not F5/SAT trials.</p><table><caption>Frozen round 3 values; no new statistics calculated by this page</caption><thead><tr><th>Curve cell</th><th>Challenger / IC reference</th></tr></thead><tbody>${Object.entries(online.per_cell).map(([cell,value])=>`<tr><td>${escape(cell)}</td><td>${f(value,6)}×</td></tr>`).join('')}<tr><td>Overall descriptive 95% interval</td><td>${f(lo,6)}–${f(hi,6)}×</td></tr></tbody></table>${link('improvement/round3/RESULTS.json','Frozen result JSON')} · <a href="https://github.com/aburan28/crypto/blob/main/docs/ic/dashboard-overview-data.json">Source hashes and overview data</a></div></details></section>

<section class="dash-section" aria-labelledby="rho-results-title"><div class="dash-section-head"><div><p class="dash-kicker">02 / Reference check</p><h2 id="rho-results-title">How does that incumbent compare with rho?</h2></div><span class="status">Historical toy panel · exploratory</span></div><div class="dash-result-grid"><figure class="dash-chart"><figcaption id="rho-graph-title"><strong>Same targets. Same online timing boundary.</strong><span>Recorded milliseconds · zero-based linear scale · lower is better</span></figcaption><div class="chart-scroll" tabindex="0" aria-label="Scrollable rho comparison graph">${rhoGraph}</div><p class="chart-note">Equal-cell geometric means of per-point three-process medians. Each arm verified 216/216 repetitions. No multi-target amortization.</p></figure><div class="dash-explanation"><span class="result-number">${f(rows.incumbent.rho_online_over_IC_online,2)}× <span>rho / IC online time</span></span><h3>A descriptive lead on this toy panel.</h3><p>Recorded IC/rho cost ratio: <strong>${f(rho.candidate_over_baseline)}×</strong>; descriptive 95% interval <strong>${f(rho.ci95[0])}–${f(rho.ci95[1])}</strong>.</p><p class="panel-summary">This is the pair-table incumbent; it does not rank F5 or SAT. Reusable preparation is excluded. The source has no host-level isolation receipt, so this wall-time ratio remains exploratory and is not a promoted CPU speedup.</p>${link('improvement/round3/README.md','Reference qualification and separate cold costs →')}</div></div></section>

<section class="dash-section" id="lab-readiness" aria-labelledby="readiness-title"><div class="dash-section-head"><div><p class="dash-kicker">03 / Solver readiness</p><h2 id="readiness-title">Correctness comes before a speed ranking.</h2></div></div><p class="dash-section-intro panel-summary">Both newly registered natural preparations reached rank 29/29 on 512 ordinary queries. F5 also recovered one registered public point. SAT has a disclosed target control and complete reusable logs, but no audited fresh target solve. Neither has an admitted fresh speed ranking.</p>${naturalStageChart}<div class="readiness"><table><caption>Evidence levels for this bounded autolab; “pending” never means zero cost</caption><thead><tr><th>Pipeline family</th><th>Disclosed target</th><th>Native execution evidence</th><th>New natural yield</th><th>Fresh F5/SAT comparison</th></tr></thead><tbody><tr><td><strong>Pair-table incumbent</strong><small>Registered round 3 reference</small></td><td><span class="status good">Verified panel</span></td><td><span class="status good">Accepted round</span></td><td>Different method<small>No F5/SAT yield inference</small></td><td><span class="status pending">Fresh paired run needed</span></td></tr><tr><td><strong>Matrix F5</strong><small>${link('native-f5-one-target-registration-v1/result-v1/RESULT.md','Audited one-target result')}</small></td><td><span class="status good">Scalar verified</span></td><td><span class="status good">Source-bound, audited</span><small>Registered n17 one-point solve</small></td><td><span class="status good">512 queries; rank 29/29</span><small>129 witnesses; 29 logs verified</small></td><td><span class="status pending">Pending</span></td></tr><tr><td><strong>SAT · CryptoMiniSat</strong><small>${link('native-sat-control-registration-v1/RESULT.md','Accepted native control')} · ${link('native-sat-million-registration-v1/result-v1/RESULT.md','Natural panel audit')}</small></td><td><span class="status good">Disclosed scalar verified</span></td><td><span class="status good">Source-bound, audited</span><small>Target control and separate natural preparation</small></td><td><span class="status good">512 queries; rank 29/29</span><small>145 witnesses; 29 logs verified; 28 unresolved</small></td><td><span class="status pending">Pending</span></td></tr></tbody></table></div>
<details class="dash-details"><summary id="control-title">Control diagnostics and the earlier SAT failure</summary><div class="dash-detail-body"><p>These are disclosed-input correctness controls. Historical Python preparation provenance remains explicit; neither native control ran new ordinary queries. Both independent audits run outside their producer intervals, so neither interval establishes the primary independently verified online speedup. Native F5 retains an outer stopwatch 500 ns below its phase sum; both values remain in the original record.</p><table><caption>Uncalibrated disclosed-input diagnostics; do not compare these as solver rankings</caption><thead><tr><th>Control</th><th>Retained producer interval</th><th>Attempts</th><th>Scope</th></tr></thead><tbody><tr><td>Historical F5</td><td>11.955969917 s</td><td>2 negatives + 1 witness</td><td>Native target worker, historical controller; ${link('prepared-f5-v3-control-v1/RESULT.md','original record')}</td></tr><tr><td>Native F5</td><td>22.420494542 s · phase sum</td><td>2 exact negatives + 1 witness</td><td>Source-bound audited execution; ${link('native-f5-control-registration-v1/RESULT.md','clock discrepancy and original audit')}</td></tr><tr><td>Native SAT</td><td>56.585480167 s</td><td>2 exact negatives + 1 witness</td><td>Source-bound audited execution; ${link('native-sat-control-registration-v1/RESULT.md','timing boundary and replay')}</td></tr><tr><td>Earlier SAT control</td><td>Not a completed online solve</td><td>8 inconclusive attempts</td><td>Failure retained; ${link('prepared-one-target-controls-v1/RESULT.md','original diagnosis')}</td></tr></tbody></table><p>Both current controls use 63 geometric points, 62 usable points before orbit folding, and 29 relation columns. These counts are distinct. Matched rho speedup for both controls remains <strong>unknown</strong>.</p></div></details></section>

<section class="dash-section" id="shared-rank-folded-table-20261003" aria-labelledby="shared-rank-folded-title">
<div class="dash-section-head"><div><p class="dash-kicker">n37 reusable setup / Stage diagnostic</p><h2 id="shared-rank-folded-title">A target-blind rank database now checks all 42 base logs.</h2></div><span class="status neutral">Correctness passed · no rho ratio</span></div>
<p class="dash-section-intro panel-summary">The frozen compact-orbit source base and signed-Frobenius three-summand table used 55 seeded <code>[a]G</code> probes, with no public Q, to collect 42 independent rows and 13 dependent rows. A separate general-curve-law verifier replayed every relation and checked all 42 solved column logs as full points. Two final producer runs have identical non-timing evidence.</p>
<div class="table-scroll"><table><caption>One target-independent n37 setup, counted group-addition equivalents. Field, hash, allocation and full memory costs remain unpriced; target descent and rho were not measured in this gate.</caption><thead><tr><th>Base</th><th>Folded table</th><th>Rank search</th><th>Elimination</th><th>Checks</th><th>Shared total</th></tr></thead><tbody><tr><td>4,260</td><td>178,710</td><td>8,498</td><td>1,110</td><td>2,248</td><td><strong>194,826</strong></td></tr></tbody></table></div>
<p class="dash-footnote">The table built 64,467 retained entries with 66,822 additions and accounts for 91.7% of the charged setup. A separately frozen follow-up recovered and independently replayed 16/16 new point-only Q directly, at a 197,454 counted cold GAE lower bound for the setup plus all 16 targets. Native work and memory remain unpriced. The separate one-target ecbench gate below pairs a fresh Q with strong rho under a new calibration; its counts are not directly subtracted from these historical counts. <a href="https://github.com/aburan28/crypto/blob/main/research/notes/ecc2k130/shared_rank_folded_table_20261003/RESULT.md">Rank evidence</a> · <a href="https://github.com/aburan28/crypto/blob/main/research/notes/ecc2k130/n37_shared_rank_targets_20261003/RESULT.md">Point-only target evidence →</a></p>
</section>

<section class="dash-section" id="n37-shared-rank-ecbench-20261004" aria-labelledby="shared-ecbench-title">
<div class="dash-section-head"><div><p class="dash-kicker">n37 one public point / paired ecbench</p><h2 id="shared-ecbench-title">The shared rank solves Q; the speed claim stays open.</h2></div><span class="status neutral">Accounting · L0 diagnostic</span></div>
<p class="dash-section-intro panel-summary">A fresh orbit-disjoint Q was solved by the 3,108-point, 42-column shared-rank IC and by single-target strong signed-Frobenius rho in all five measured rounds. An independent Linux x86-64 auditor reproduced all 15 measured records exactly. The target-independent table and 42 verified base logs were built before each IC online clock; each target was a direct m3 hit. Both arms have unpriced native work, and the Mac wall times remain L0 diagnostics.</p>
<div class="table-scroll"><table><caption>Same Q, five verified one-target repetitions. Every S and ratio is a lower-bound or bounded cold diagnostic; online times are L0 observations, not an accepted speedup.</caption><thead><tr><th>Variant</th><th>Verified</th><th>Mean cold S</th><th>S / floor</th><th>S / rho</th><th>Median target online</th></tr></thead><tbody><tr><td>Strong signed-Frobenius rho</td><td>5/5</td><td>${f(sharedComparison.a.mean_s,3)}</td><td>${f(sharedComparison.a.mean_ratio_to_floor,3)}</td><td>1.000</td><td>${f(sharedEcbench.rho_online_ns.median_ns/1000,3)} µs</td></tr><tr><td>Shared-rank IC</td><td>5/5</td><td>${f(sharedComparison.b.mean_s,3)}</td><td>${f(sharedComparison.b.mean_ratio_to_floor,3)}</td><td>${f(sharedComparison.ops.ratio_b_over_a,3)}</td><td>${f(sharedEcbench.ic_online_ns.median_ns/1000,3)} µs</td></tr><tr><td>Shared-rank IC A/A control</td><td>5/5</td><td>${f(sharedControl.b.mean_s,3)}</td><td>${f(sharedControl.b.mean_ratio_to_floor,3)}</td><td>${f(sharedComparison.ops.ratio_b_over_a,3)}</td><td>${f(sharedEcbench.control_online_ns.median_ns/1000,3)} µs</td></tr></tbody></table></div>
<p class="dash-footnote">Five repeats on one Q quantify run noise, not target variation. The A/A IC online ratio ranged 0.833–1.287. An independent Linux x86-64 auditor reproduced ${sharedIndependent.replays.length}/${sharedIndependent.replays.length} measured records on another host class; this checks deterministic answers and counts, not the L0 wall times. The eight-point column sweep below addresses target variation in the same n37 method, but still needs isolated timing and common native-work pricing. No n131 transfer follows. <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_shared_rank_20261004/independent_validation_20261004/RESULT.md">Independent replay</a> · <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_shared_rank_20261004/README.md">Protocol, session, audit and decision</a> · <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_shared_rank_20261004/RESULT_FINAL.json">Native analysis →</a></p>
</section>

<section class="dash-section" id="n37-rank-columns-20261004" aria-labelledby="rank-columns-title">
<div class="dash-section-head"><div><p class="dash-kicker">n37 / eight one-target workloads / factor-base sweep</p><h2 id="rank-columns-title">16 folded columns cut counted cold work; rho remains ahead in this diagnostic.</h2></div><span class="status neutral">Independent replay · L0 counts</span></div>
<p class="dash-section-intro panel-summary">The preregistered K={4,8,12,16,24,42} sweep solved the same eight public points with each base and strong signed-Frobenius rho, five times per point. All 320 measured executions verified and replayed exactly on another host class. K16 has 1,184 actual usable points and rank 16/16. Its cold counted GAE is 0.349 of K42, 95% interval [0.343, 0.359], passing the frozen engineering gate; K8 was statistically close in that incomplete unit. A separate whole-solve Callgrind census confirms K16/K42 at ${f(callgrind.comparisons[0].ratio_of_sums)} [${f(callgrind.comparisons[0].bootstrap_95[0])}, ${f(callgrind.comparisons[0].bootstrap_95[1])}] but finds K16/rho at ${f(callgrind.comparisons[1].ratio_of_sums)} [${f(callgrind.comparisons[1].bootstrap_95[0])}, ${f(callgrind.comparisons[1].bootstrap_95[1])}] in that separate instruction unit. The untouched K8/K16 confirmation follows below; isolated online wall speed remains unknown.</p>
<div class="table-scroll"><table><caption>Cold counted GAE lower-bound S per √r; 40 verified one-target runs per arm on the same Q panel. The IC/rho quotient is a ratio of incomplete cost estimates and does not bound true speed.</caption><thead><tr><th>Arm</th><th>Usable points</th><th>Folded columns</th><th>Mean cold S</th><th>Counted IC / rho diagnostic</th></tr></thead><tbody>
<tr><td>Strong signed-Frobenius rho</td><td>—</td><td>—</td><td>${f(rankColumns.reference.mean_cold_s_lower_bound,3)}</td><td>1.000</td></tr>
${rankColumns.candidates.map(row => `<tr><td>${escape(row.arm)}</td><td>${row.usable_points.toLocaleString('en-US')}</td><td>${row.folded_columns}</td><td>${f(row.mean_cold_s_lower_bound,3)}</td><td>${f(row.cold_counted_over_rho_diagnostic,3)}</td></tr>`).join('\n')}
</tbody></table></div>
<p class="dash-footnote">K16's counted preparation mean is 31,312.49 GAE and its target-dependent mean is 748.32 GAE; K42 is 91,712.76 plus 156.89. The earlier K16/K8 counted ratio was 0.990 [0.927, 1.057]; the fresh confirmation below uses new targets and prices native implementation work. The identical K42 control had a counted A/A ratio of 1.000. The Mac producer is L0 and field arithmetic, hashing, allocation and modular combination are unpriced in both methods; online wall speedup and fully priced cold speedup remain unset. The independent Linux receipt reproduces ${rankColumnsAudit.replays.length}/${rankColumnsAudit.replays.length} measured records, not an isolated producer timing. <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_rank_columns_20261004/RESULT.md">Decision and limitations</a> · <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_rank_columns_20261004/RESULT.json">Native analysis</a> · <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_rank_columns_20261004/independent_validation/README.md">Replay evidence →</a></p>
<p class="dash-footnote">A follow-on whole-solve Callgrind census on the same eight public points counts field, hash, allocation and other user-space work inside the implementation: K16/K42 is ${f(callgrind.comparisons[0].ratio_of_sums)} [${f(callgrind.comparisons[0].bootstrap_95[0])}, ${f(callgrind.comparisons[0].bootstrap_95[1])}], but K16 still costs ${f(callgrind.comparisons[1].ratio_of_sums)} [${f(callgrind.comparisons[1].bootstrap_95[0])}, ${f(callgrind.comparisons[1].bootstrap_95[1])}] times strong rho in this separate instruction unit. All ${callgrind.profiles} profiles reproduced the archived logs and the K42 A/A maximum difference was ${(100*callgrind.aa_max_absolute_relative_deviation).toFixed(4)}%. This is neither native wall time nor the primary online speed metric. <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_callgrind_solve_20261004/RESULT.md">Frozen instruction decision and raw replay →</a></p>
</section>

<section class="dash-section" id="n37-k8-k16-20261004" aria-labelledby="k8-confirm-title">
<div class="dash-section-head"><div><p class="dash-kicker">n37 / 16 untouched one-target workloads / complete solve instructions</p><h2 id="k8-confirm-title">K8 wins cold implementation cost; target-only attribution follows.</h2></div><span class="status neutral">Accounting · independent replay · L0</span></div>
<p class="dash-section-intro panel-summary">The fresh public-point panel verified 384/384 native executions and 64/64 Callgrind profiles, with all 320 measured native rows independently replayed on Linux. K16/K8 whole-solve instructions are ${f(k8Confirm.comparisons[0].ratio_of_sums)} [${f(k8Confirm.comparisons[0].bootstrap_95[0])}, ${f(k8Confirm.comparisons[0].bootstrap_95[1])}], clearing the preregistered 5% K8 cold-route gate. K8/rho remains ${f(k8Confirm.comparisons[1].ratio_of_sums)} [${f(k8Confirm.comparisons[1].bootstrap_95[0])}, ${f(k8Confirm.comparisons[1].bootstrap_95[1])}] in the same simulated instruction unit. The new target-only instruction result appears below; the primary isolated online wall result remains unset.</p>
<div class="table-scroll"><table><caption>Complete ecbench method-solve Callgrind Ir on the same 16 public points; one profile per arm and target. The cold method-solve boundary includes reusable setup and is not the online interval.</caption><thead><tr><th>Arm</th><th>Usable points</th><th>Folded columns</th><th>Verified</th><th>Mean solve Ir</th><th>Ir / √r</th></tr></thead><tbody>
<tr><td>Strong signed-Frobenius rho</td><td>—</td><td>—</td><td>16/16</td><td>${Math.round(k8Confirm.arms.find(row => row.arm === 'rho-strong').mean_ir).toLocaleString('en-US')}</td><td>${f(k8Confirm.arms.find(row => row.arm === 'rho-strong').s_ir_per_target,2)}</td></tr>
<tr><td>K8</td><td>592</td><td>8</td><td>16/16</td><td>${Math.round(k8Confirm.arms.find(row => row.arm === 'ic-k8').mean_ir).toLocaleString('en-US')}</td><td>${f(k8Confirm.arms.find(row => row.arm === 'ic-k8').s_ir_per_target,2)}</td></tr>
<tr><td>K16</td><td>1,184</td><td>16</td><td>16/16</td><td>${Math.round(k8Confirm.arms.find(row => row.arm === 'ic-k16').mean_ir).toLocaleString('en-US')}</td><td>${f(k8Confirm.arms.find(row => row.arm === 'ic-k16').s_ir_per_target,2)}</td></tr>
<tr><td>Identical K16 control</td><td>1,184</td><td>16</td><td>16/16</td><td>${Math.round(k8Confirm.arms.find(row => row.arm === 'ic-k16-control').mean_ir).toLocaleString('en-US')}</td><td>${f(k8Confirm.arms.find(row => row.arm === 'ic-k16-control').s_ir_per_target,2)}</td></tr>
</tbody></table></div>
<p class="dash-footnote">The older incomplete counted-cost estimate had K16/K8 = 1.0175 [0.973, 1.059] on these new targets and could not choose a base. K8's pair table has 2,344 entries versus K16's 9,412, while K16 averaged 477 counted target GAE versus K8's 3,385. The largest paired K16 A/A deviation was ${(100*k8Confirm.aa_max_absolute_relative_deviation).toFixed(4)}%. These phase counters explain why both bases must remain in the n41/n53 online study: the cold instruction choice cannot settle target-only cost. <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_k8_k16_20261004/RESULT.md">Decision, raw profile archive, exact IC1 claims and limits →</a></p>
</section>

<section class="dash-section" id="n37-online-ir-20261004" aria-labelledby="n37-online-ir-title">
<div class="dash-section-head"><div><p class="dash-kicker">n37 / 16 new public one-target workloads / target-only instructions</p><h2 id="n37-online-ir-title">K16 cuts target work; an isolated wall test is next.</h2></div><span class="status neutral">Accounting · independent replay · L0</span></div>
<p class="dash-section-intro panel-summary">The fresh 16-point panel verified ${onlineIr.profiles}/${onlineIr.profiles} Callgrind profiles and independently replayed ${onlineIrReceipt.replays.length}/${onlineIrReceipt.replays.length} measured native executions. K8/K16 target-only instructions are ${f(onlineIr.online.comparisons[0].ratio_of_sums)} [${f(onlineIr.online.comparisons[0].bootstrap_95[0])}, ${f(onlineIr.online.comparisons[0].bootstrap_95[1])}], above the frozen 1.10 K16-priority gate. K16/rho target-only instructions are ${f(onlineIr.online.comparisons[2].ratio_of_sums)} [${f(onlineIr.online.comparisons[2].bootstrap_95[0])}, ${f(onlineIr.online.comparisons[2].bootstrap_95[1])}]. This is simulated instruction attribution, not an isolated online wall speedup.</p>
<div class="table-scroll"><table><caption>One-target Callgrind Ir from the frozen solve and target interval markers; 16 same-point verified profiles per arm.</caption><thead><tr><th>Arm</th><th>Actual usable points</th><th>Folded columns</th><th>Mean target-only Ir</th><th>Mean complete-solve Ir</th></tr></thead><tbody>
${[['rho-strong','Strong signed-Frobenius rho','—','—'],['ic-k8','K8','592','8'],['ic-k16','K16','1,184','16'],['ic-k16-control','Identical K16 control','1,184','16']].map(([arm,label,points,columns])=>`<tr><td>${label}</td><td>${points}</td><td>${columns}</td><td>${Math.round(onlineIr.online.arms.find(row=>row.arm===arm).mean_ir).toLocaleString('en-US')}</td><td>${Math.round(onlineIr.complete_solve.arms.find(row=>row.arm===arm).mean_ir).toLocaleString('en-US')}</td></tr>`).join('\n')}
</tbody></table></div>
<p class="dash-footnote">The cold choice reverses: K8/K16 complete-solve Ir is ${f(onlineIr.complete_solve.comparisons[0].ratio_of_sums)} [${f(onlineIr.complete_solve.comparisons[0].bootstrap_95[0])}, ${f(onlineIr.complete_solve.comparisons[0].bootstrap_95[1])}], as K16's larger base and pair table cost more reusable work. The target-only K16 A/A difference is exactly zero. Both bases stay in the n41/n53 study; no ECC2K-130 transfer follows from n37. The Mac wall session earned L0, so the primary online speedup stays unknown. <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_online_ir_20261004/RESULT.md">Decision, raw archive, exact IC1 claims, and limits →</a></p>
</section>

<section class="dash-section" id="n37-native-online-wall-20261004" aria-labelledby="n37-native-online-wall-title">
<div class="dash-section-head"><div><p class="dash-kicker">n37 / 16 new public one-target workloads / native wall screen</p><h2 id="n37-native-online-wall-title">K16's native online lead remains exploratory.</h2></div><span class="status neutral">L1 · 320/320 independent replays · noise gate failed</span></div>
<p class="dash-section-intro panel-summary">All 384 executions recovered and verified their public-point logarithms. The frozen primary target, <code>${escape(nativeWall.primary_one_target.workload_id)}</code>, has descriptive rho/K16 online wall ${f(nativeWall.primary_one_target.descriptive_rho_over_k16)}×, but no admitted speedup: all measured rows earned L1, and the identical K16 A/A maximum deviation was ${f(100*nativeWall.k16_aa_max_relative_deviation,2)}%, above the frozen 5% gate.</p>
<div class="table-scroll"><table><caption>Primary target-zero mean of five paired one-target rounds; hosted wall milliseconds are exploratory.</caption><thead><tr><th>Arm</th><th>Usable points</th><th>Folded columns</th><th>Online wall, ms</th><th>Cold solve wall, ms (16-target panel)</th></tr></thead><tbody>
${[['rho-strong','Strong signed-Frobenius rho','—','—'],['ic-k8','K8','592','8'],['ic-k16','K16','1,184','16'],['ic-k16-control','Identical K16 control','1,184','16']].map(([arm,label,points,columns])=>`<tr><td>${label}</td><td>${points}</td><td>${columns}</td><td>${f(nativeWall.primary_one_target.online_ns_totals[arm]/5e6,6)}</td><td>${f(nativeWall.arm_stage_totals[arm].cold_solve_ns_total/80e6,6)}</td></tr>`).join('\n')}
</tbody></table></div>
<p class="dash-footnote">Across the 16 separate one-target workloads, the secondary descriptive K8/K16 online ratio is ${f(nativeWall.panel_descriptive_ratios.k8_over_k16.ratio_of_sums)} [${f(nativeWall.panel_descriptive_ratios.k8_over_k16.target_block_bootstrap_95[0])}, ${f(nativeWall.panel_descriptive_ratios.k8_over_k16.target_block_bootstrap_95[1])}]. The frozen decision is <code>${escape(nativeWall.decision)}</code>: carry both bases to an auditable L2 host; retain K8 and K16 at n41/n53. The page's primary rho/IC online speedup remains unknown, and this n37 screen supplies no ECC2K-130 transfer. <a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_n37_native_online_wall_20261004/RESULT.md">Raw records, independent replay, exact IC1 identities and decision →</a></p>
</section>

<section class="dash-section" id="lab-pipeline" aria-labelledby="pipeline-title"><div class="dash-section-head"><div><p class="dash-kicker">04 / Pipeline map</p><h2 id="pipeline-title">What the solver actually has to do.</h2></div></div><p class="dash-section-intro panel-summary">Reusable preparation produces factor-base logs. The online solve consumes one registered public point and includes every target-dependent attempt through independent scalar verification.</p>
<div class="pipeline-zone"><div class="pipeline-zone-label">Reusable preparation <span>Cost and memory reported separately</span></div><ol class="pipeline"><li><b>1</b><strong>Build the factor base</strong><span>Choose exact subgroup points and declare sign / Frobenius orbit folding.</span></li><li><b>2</b><strong>Find and verify relations</strong><span>Ordinary queries → point decomposition → verified, useful matrix rows.</span></li><li><b>3</b><strong>Solve the relation matrix</strong><span>Check rank and recover verified factor-base logarithms.</span></li></ol></div>
<div class="pipeline-zone online-zone"><div class="pipeline-zone-label">One target · primary online clock <span>Start at target-dependent work → stop after independent verification</span></div><ol class="pipeline"><li><b>4</b><strong>Decompose the target</strong><span>Generate queries, encode and solve. Charge failed and timed-out attempts.</span></li><li><b>5</b><strong>Recover its logarithm</strong><span>Check the relation, perform target descent and combine the known factor logs.</span></li><li><b>6</b><strong>Verify the answer</strong><span>Independently check that the recovered scalar maps the generator to this target.</span></li></ol></div>
<div class="stage-details"><details><summary>Factor bases &amp; Frobenius orbits</summary><p>Record the actual usable point set before folding, then folded columns and final rank separately. Different bases are different complete candidates even when their counts match.</p></details><details><summary>PDP: F4 / F5 / SAT</summary><p>PDP means point decomposition. Its internal polynomial or Macaulay matrix belongs to this stage. Solvers must return a verified point relation, not merely a satisfiable encoding. Every inconclusive attempt stays visible.</p></details><details><summary>Final linear algebra &amp; descent</summary><p>The final relation matrix is distinct from the solver’s internal matrix. Rank, modular arithmetic, recovered logs, recursive target work and the final scalar check all need evidence.</p></details></div><p class="dash-footnote">Rho must solve the same one point with the same resource limits. Fixture generation, process launch and input loading are outside both online intervals. A missing phase cost remains unknown.</p></section>

<section class="dash-section" id="lab-next" aria-labelledby="next-title"><div class="dash-section-head"><div><p class="dash-kicker">05 / Tournament design</p><h2 id="next-title">A local win is a candidate, not a final answer.</h2></div></div><div class="tournament-flow" aria-label="Proposed tournament flow"><article><span class="flow-number">01 / EXPLORE</span><strong>Keep diverse candidates</strong><p>Different bases, solver families and stage combinations. Preserve alternatives with different resource tradeoffs.</p></article><article><span class="flow-number">02 / ADMIT</span><strong>Pass correctness gates</strong><p>Source identity, natural yield, failed attempts, complete costs and verified recovered answers.</p></article><article><span class="flow-number">03 / COMPARE</span><strong>Pair frozen workloads</strong><p>Same fresh point, resources and reference. Screening and untouched confirmation have separate roles.</p></article><article><span class="flow-number">04 / DECIDE</span><strong>Promote or retain</strong><p>Apply the declared uncertainty gate. Keep negative outcomes and useful alternatives for later combinations.</p></article></div><div class="next-callout"><span class="status pending">Next gate</span><div><strong>Validate prepared CMS transport, then freeze the four-arm target.</strong><p>The SAT natural log table is complete. The marked CMS binary still needs its predeclared disclosed-input parity audit; F5, SAT, incumbent and rho then need prepublished source, one fresh point and one-use runs. The three historical confirmation sets remain closed.</p></div></div><p class="dash-footnote panel-summary">This diagram describes the proposed process. It does not claim that a fresh F5/SAT tournament has run, or that any implementation is globally fastest.</p></section>

<section class="dash-library-intro" id="evidence-search" aria-labelledby="library-title"><p class="dash-kicker">06 / Evidence</p><h2 id="library-title">Open the detail you need.</h2><p class="panel-summary">The historical ledger contains other curves, cost units and questions. Its reports are preserved; they are not a combined leaderboard.</p><label for="evidence-query">Search historical reports</label><input id="evidence-query" type="search" placeholder="Frobenius, F5, rho, factor base…" autocomplete="off" aria-controls="evidence-results"><p id="evidence-count" aria-live="polite"></p><ul id="evidence-results"></ul><noscript><p>Open the full ledger below and use your browser’s Find command.</p></noscript></section><details class="dash-details" id="historical-regimes"><summary>Historical best results by regime · different workloads and accounting</summary><div class="dash-detail-body">${historical}</div></details>
${historicalDetails}
</main>
`;
page = replace(page, START, END, front);
page = replace(page, CSS_START, CSS_END, '\n' + readFileSync(resolve(ROOT,'docs/ic/dashboard-overview.css'),'utf8') + '\n');
page = replace(page, '<script id="ic-overview-interactions">', '</script>', '\n' + readFileSync(resolve(ROOT,'docs/ic/dashboard-overview-interactions.js'),'utf8') + '\n');
page = page.replace(/<title>[^<]*<\/title>/, '<title>IC Research Lab — Decisions, Pipeline &amp; Evidence</title>');
page = page.replace(/<meta name="description" content="[^"]*">/, '<meta name="description" content="Understand the IC tournament decision, solver readiness, one-target pipeline and evidence. Frozen toy-panel measurements, visible uncertainty and preserved failures.">');
// The self-contained dashboard should work offline without a font request.
page = page.replace(/<link[^>]*href="https:\/\/fonts\.[^"]*"[^>]*>\n?/g, '');
page = page.replace('it loads two webfont\n  families from Google Fonts and needs nothing else.', 'its overview needs no network requests\n  or installed dependencies.');
assert.equal(between(page, '<div id="legacy-evidence">', '<!-- END IC EVIDENCE LIBRARY -->'), ledgerBefore, 'Historical ledger changed');
assert.ok(page.includes(historical), 'Historical regime summary lost');
for (const panel of retainedPanels) assert.ok(page.includes(panel.html), `Historical panel lost: ${panel.id}`);
const data = {
  schema_version: 3, scope: 'bounded IC autolab overview; not all repository research',
  sources: [roundPin,f5Pin,nativeF5Pin,oldSatPin,satPin,ordinaryPin,satMillionPin,f5TargetPin,pilotNegativePin,pilotFeasiblePin,sharedRankPin,solverComparisonPin,sharedTargetsPin,sharedEcbenchPin,sharedControlPin,sharedIndependentPin,sharedIndependentClaimPin,sharedIndependentCheckPin,rankColumnsPin,rankColumnsAuditPin,callgrindPin,k8ConfirmPin,onlineIrPin,onlineIrReceiptPin,nativeWallPin,nativeWallEvidencePin], confirmation_online: online,
  confirmation_familywise_online_upper: upper, paired_targets: confirmation.paired_cases,
  historical_f5_control: f5, native_f5_control: nativeF5, native_sat_control: sat, historical_incomplete_sat_control: oldSat.arms.sat,
  native_ordinary_stage: {ordinary_queries:ordinary.ordinary_queries,
    f5:{outcome_mix:ordinary.f5.outcome_mix,final_rank:ordinary.f5.final_rank,folded_columns:ordinary.f5.folded_columns,independently_recovered_logs:ordinary.f5.independently_recovered_logs},
    cms:{outcome_mix:ordinary.cms.outcome_mix,final_rank:ordinary.cms.final_rank,folded_columns:ordinary.cms.folded_columns,feasible_inconclusive_queries:ordinary.cms.feasible_inconclusive_queries},
    online_speedup:ordinary.online_speedup},
  native_sat_million_stage: {ordinary_queries:satMillion.audited_ordinary_queries,
    outcome_mix:satMillion.mathematics.ordinary_outcome_mix,
    accepted_rows:satMillion.mathematics.accepted_rows,
    final_rank:satMillion.mathematics.rank,
    folded_columns:satMillion.mathematics.folded_columns,
    independently_recovered_logs:satMillion.mathematics.independently_recovered_column_logs.length,
    source_bound_preparation_admitted:satMillion.source_bound_execution_admitted,
    target_solve_admitted:false,online_speedup:null},
  native_f5_one_target: {source_bound_execution_admitted:f5Target.source_bound_execution_admitted,
    verified_recovery:f5Target.verified_recovery,recorded_online_interval_ns:f5Target.recorded_online_interval_ns,
    headline_eligible:f5Target.headline_eligible,online_speedup:f5Target.online_speedup},
  disclosed_sat_budget_pilot: {negative_status:pilotNegative.cold.status,
    feasible_status:pilotFeasible.cold.status,source_bound_target_admitted:false,
    natural_yield_estimate:null},
  rho_online_comparison: rho,
  shared_rank_ecbench: sharedEcbench,
  shared_rank_independent_replay: {session_id:sharedIndependent.session_id, auditor_env_class_id:sharedIndependent.auditor_env_class_id, replays:sharedIndependent.replays.length, receipt_sha256:sharedIndependentSha, claim_check_status:sharedIndependentCheck.status},
  n37_rank_columns: rankColumns,
  n37_callgrind: callgrind,
  n37_k8_k16_callgrind: k8Confirm,
  n37_online_ir: onlineIr,
  n37_online_ir_independent_replay: {session_id:onlineIrReceipt.session_id, auditor_env_class_id:onlineIrReceipt.auditor_env_class_id, replays:onlineIrReceipt.replays.length, receipt_sha256:onlineIrReceiptSha},
  n37_native_online_wall: nativeWall,
  n37_native_online_wall_evidence: nativeWallEvidence,
  n37_rank_columns_independent_replay: {session_id:rankColumnsAudit.session_id, auditor_env_class_id:rankColumnsAudit.auditor_env_class_id, replays:rankColumnsAudit.replays.length, receipt_sha256:rankColumnsAuditSha},
  rho_online_table: Object.values(rows).map(row=>Object.fromEntries(['alias','online_ms','rho_online_over_IC_online','verified','scheduled'].map(key=>[key,row[key]]))),
  historical_ledger_sha256: sha(ledgerBefore), historical_regime_summary_sha256: sha(historical),
  historical_overview_panels: retainedPanels.map(panel=>({section_id:panel.id,sha256:sha(panel.html)})),
  progress_timeline_sha256: sha(progressRaw)
};
const output = JSON.stringify(data,null,2) + '\n';
if (process.argv.includes('--check')) {
  assert.equal(readFileSync(PAGE,'utf8'), page, 'Dashboard is stale; run node tools/render_ic_dashboard.mjs');
  assert.equal(readFileSync(resolve(ROOT,'docs/ic/dashboard-overview-data.json'),'utf8'), output, 'Dashboard data is stale');
  console.log('PASS: evidence pins, status gates, generated page/data and unchanged historical ledger');
} else {
  writeFileSync(PAGE,page);
  writeFileSync(resolve(ROOT,'docs/ic/dashboard-overview-data.json'),output);
  console.log('Rendered source-pinned dashboard; historical ledger and regime summary preserved.');
}
