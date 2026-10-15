// Compare immutable repository objects; this does not execute research code.
import {execFileSync} from 'node:child_process';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
const git = (...args) => execFileSync('git', args, {encoding: 'utf8', maxBuffer: 32 * 1024 * 1024}).trim();
const prior = '6351a6a3df065c5ffc7b81aba5d4a82157b022e6';
const incoming = 'f2c9202ca48ed16f086815fda9c0ab73f2aabcda';
const tree = git('write-tree');
function entries(ref) {
  return new Map(git('ls-tree', '-r', ref).split('\n').map(line => {
    const [record, path] = line.split('\t');
    return [path, record];
  }));
}
const merged = entries(tree);
const studies = ['research/ecbench_open_questions_20261014/', 'research/ecbench_leads_20261015/'];
const ownFrozen = path => studies.some(prefix => path.startsWith(prefix)) ||
  path.startsWith('benchmarks/20261014-ecbench-ten-questions/') ||
  path.startsWith('benchmarks/20261015-ecbench-leads/');
const incomingFrozen = path => (path.startsWith('research/') && !studies.some(prefix => path.startsWith(prefix))) ||
  (path.startsWith('benchmarks/') && path !== 'benchmarks/README.md' && !ownFrozen(path)) ||
  path.startsWith('docs/ic/factor-base-encyclopedia/') || path.startsWith('reports/') || path.startsWith('research_notes/');
const results = [];
for (const [label, ref, select] of [['prior-campaign', prior, ownFrozen], ['incoming-evidence', incoming, incomingFrozen]]) {
  let checked = 0;
  for (const [path, object] of entries(ref)) {
    if (!select(path)) continue;
    if (merged.get(path) !== object) throw new Error(`${label}: changed or absent ${path}`);
    checked++;
  }
  results.push({label, ref, checked, identical: true});
}
// Every non-index incoming library file is exactly its main-parent object.
for (const path of ['src/cryptanalysis/ecbench/fp4.rs', 'src/cryptanalysis/ecbench/workload.rs', 'src/cryptanalysis/koblitz_groebner.rs', 'src/cryptanalysis/koblitz_index_calculus.rs']) {
  if (merged.get(path) !== entries(incoming).get(path)) throw new Error(`incoming source changed: ${path}`);
}
const index = readFileSync('benchmarks/README.md', 'utf8').split('\n').filter(line => line.startsWith('|'));
for (const ref of [prior, incoming]) {
  for (const row of git('show', `${ref}:benchmarks/README.md`).split('\n').filter(line => line.startsWith('|'))) {
    if (!index.includes(row)) throw new Error(`benchmark row lost from ${ref}`);
  }
}
const measured = createHash('sha256').update(readFileSync('/workspace/measured-bin/ecbench-leads-0cbc')).digest('hex');
if (measured !== '50a5739fadcb08ab8f170d4fe03ba3427fb0fd387c1a5918843fcda4d68ff1fb') throw new Error('measured binary changed');
console.log(JSON.stringify({tree_before_validation_receipts: tree, comparisons: results, benchmark_rows_from_both_parents_preserved: true, measured_binary_sha256: measured, scientific_campaign_executed: false}, null, 2));
