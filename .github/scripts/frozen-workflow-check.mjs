#!/usr/bin/env node
// Historical replay manifests bind workflow metadata as well as study inputs.
// Scheduling may change without rewriting that record: execute the unchanged
// validator in a copied current source view with only its pinned workflow bytes.
import { execFileSync, spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { copyFileSync, existsSync, lstatSync, mkdirSync, mkdtempSync, readFileSync, readdirSync, rmSync, symlinkSync, writeFileSync } from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const sha = bytes => createHash('sha256').update(bytes).digest('hex');

function relative(file) {
  if (typeof file !== 'string' || !file || path.posix.isAbsolute(file)
      || file.includes('\\') || file.split('/').includes('..') || file.split('/').includes('.git')) {
    throw new Error(`unsafe source-view path: ${file}`);
  }
  return path.posix.normalize(file);
}

function copyCurrent(cwd, view, file) {
  const source = path.join(cwd, file);
  const destination = path.join(view, file);
  const info = lstatSync(source); // Missing current inputs must fail, never roll back.
  if (info.isSymbolicLink()) throw new Error(`symbolic source-view input: ${file}`);
  if (info.isDirectory()) {
    mkdirSync(destination, { recursive: true });
    for (const child of readdirSync(source)) copyCurrent(cwd, view, path.posix.join(file, child));
  } else if (info.isFile()) {
    mkdirSync(path.dirname(destination), { recursive: true });
    copyFileSync(source, destination);
  } else {
    throw new Error(`unsupported source-view input: ${file}`);
  }
}

function regularBytes(root, file, label) {
  const target = path.join(root, file);
  const info = lstatSync(target);
  if (!info.isFile() || info.isSymbolicLink()) throw new Error(`${label} is not a regular file: ${file}`);
  return readFileSync(target);
}

function ensureCommit(cwd, revision) {
  if (!/^[0-9a-f]{40,64}$/.test(revision ?? '')) throw new Error('expected an immutable Git revision');
  try {
    execFileSync('git', ['cat-file', '-e', `${revision}^{commit}`], { cwd, stdio: ['ignore', 'pipe', 'pipe'] });
  } catch {
    execFileSync('git', ['fetch', '--filter=blob:none', '--no-tags', '--depth=1', 'origin', revision],
      { cwd, stdio: ['ignore', 'pipe', 'pipe'] });
  }
}

function applyPinnedInputView({ cwd, view, manifest }) {
  manifest = relative(manifest);
  const configBytes = regularBytes(cwd, manifest, 'pinned-input manifest');
  const config = JSON.parse(configBytes.toString('utf8'));
  if (config.schema !== 'frozen-input-view/v1' || !Array.isArray(config.inputs) || !config.inputs.length) {
    throw new Error('invalid pinned-input manifest schema');
  }
  const overlays = [];
  for (const item of config.inputs) {
    const target = relative(item.target);
    const revision = item.revision;
    const sourcePath = relative(item.source_path);
    const originalManifest = relative(item.original_manifest);
    const sourceKey = item.original_source_key;
    const expected = item.sha256;
    if (typeof sourceKey !== 'string' || !/^[0-9a-f]{64}$/.test(expected ?? '')
        || !/^[0-9a-f]{64}$/.test(item.original_manifest_sha256 ?? '')) {
      throw new Error(`invalid pinned-input digests: ${target}`);
    }
    const originalBytes = regularBytes(cwd, originalManifest, 'original freeze manifest');
    if (sha(originalBytes) !== item.original_manifest_sha256) {
      throw new Error(`original freeze manifest changed: ${originalManifest}`);
    }
    const viewOriginalBytes = regularBytes(view, originalManifest, 'view original freeze manifest');
    if (!viewOriginalBytes.equals(originalBytes)) throw new Error(`view changed original freeze manifest: ${originalManifest}`);
    const frozen = JSON.parse(originalBytes.toString('utf8'));
    if (frozen.source_sha256?.[sourceKey] !== expected) {
      throw new Error(`pinned input does not match original frozen source: ${sourceKey}`);
    }
    const live = regularBytes(cwd, target, 'live pinned-input target');
    const currentSources = {};
    if (!Array.isArray(item.current_equal) || !item.current_equal.some(pair =>
      Array.isArray(pair) && pair.length === 2 && pair.includes(target))) {
      throw new Error(`pinned input lacks a current-source equality check: ${target}`);
    }
    for (const pair of item.current_equal) {
      if (!Array.isArray(pair) || pair.length !== 2) throw new Error(`invalid current-source pair: ${target}`);
      const [left, right] = pair.map(relative);
      if (left === right || pair.filter(name => relative(name) === target).length !== 1) {
        throw new Error(`invalid current-source equality pair: ${target}`);
      }
      const leftBytes = regularBytes(cwd, left, 'current source');
      const rightBytes = regularBytes(cwd, right, 'current source');
      currentSources[left] = sha(leftBytes);
      currentSources[right] = sha(rightBytes);
      if (!leftBytes.equals(rightBytes)) throw new Error(`current source inputs differ: ${left} != ${right}`);
    }
    ensureCommit(cwd, revision);
    const historical = execFileSync('git', ['show', `${revision}:${sourcePath}`],
      { cwd, maxBuffer: 16 * 1024 * 1024, stdio: ['ignore', 'pipe', 'pipe'] });
    if (sha(historical) !== expected) throw new Error(`pinned Git source hash mismatch: ${sourcePath}`);
    const viewTarget = path.join(view, target);
    if (!existsSync(viewTarget) || !regularBytes(view, target, 'view pinned-input target').equals(live)) {
      throw new Error(`current source view missing or changed before overlay: ${target}`);
    }
    writeFileSync(viewTarget, historical);
    if (!regularBytes(view, target, 'overlaid pinned-input target').equals(historical)) {
      throw new Error(`historical input overlay failed: ${target}`);
    }
    overlays.push({ target, revision, source_path: sourcePath, sha256: expected,
      original_manifest: originalManifest, original_manifest_sha256: item.original_manifest_sha256,
      original_source_key: sourceKey, live_target_sha256: sha(live), current_sources: currentSources,
      view_sha256: sha(historical) });
  }
  return overlays;
}

export function withFrozenWorkflow({ cwd = process.cwd(), workflow, revision, manifest, inputViewManifest,
  roots = [], command }) {
  workflow = relative(workflow);
  manifest = relative(manifest);
  roots = roots.map(relative);
  ensureCommit(cwd, revision);
  if (!Array.isArray(command) || !command.length || command.some(arg => typeof arg !== 'string')) {
    throw new Error('expected a command and argument array');
  }
  const fromGit = file => execFileSync('git', ['show', `${revision}:${file}`],
    { cwd, maxBuffer: 16 * 1024 * 1024, stdio: ['ignore', 'pipe', 'pipe'] });
  const manifestBytes = readFileSync(path.join(cwd, manifest));
  if (!manifestBytes.equals(fromGit(manifest))) throw new Error('original frozen manifest changed');
  const frozen = JSON.parse(manifestBytes.toString('utf8'));
  const workflowBytes = fromGit(workflow);
  if (!/^[0-9a-f]{64}$/.test(frozen.files?.[workflow] ?? '') || sha(workflowBytes) !== frozen.files[workflow]) {
    throw new Error('historical workflow does not match its frozen digest');
  }
  const liveWorkflowSha = sha(readFileSync(path.join(cwd, workflow)));
  const ownInputs = Object.entries(frozen.files).map(([name, digest]) => {
    const file = relative(name.startsWith('.github/') ? name : path.posix.join(path.posix.dirname(manifest), name));
    if (!/^[0-9a-f]{64}$/.test(digest)) throw new Error(`invalid frozen digest: ${name}`);
    return { file, digest };
  });
  const view = mkdtempSync(path.join(os.tmpdir(), 'crypto-frozen-workflow-'));
  try {
    // Real copies prevent validators/tests from altering the live checkout.
    for (const root of roots) copyCurrent(cwd, view, root);
    copyCurrent(cwd, view, manifest);
    for (const { file } of ownInputs) if (file !== workflow) copyCurrent(cwd, view, file);
    mkdirSync(path.dirname(path.join(view, workflow)), { recursive: true });
    writeFileSync(path.join(view, workflow), workflowBytes);
    for (const { file, digest } of ownInputs) {
      if (sha(readFileSync(path.join(view, file))) !== digest) throw new Error(`frozen input drift: ${file}`);
    }
    const pinnedInputViews = inputViewManifest
      ? applyPinnedInputView({ cwd, view, manifest: inputViewManifest }) : [];
    // Original preflights retain their Git ancestry assertions, against the
    // current checkout's unchanged object database; no historical code checkout.
    const gitPath = path.join(cwd, '.git');
    if (!existsSync(gitPath)) throw new Error('current Git metadata is missing');
    symlinkSync(gitPath, path.join(view, '.git'));
    console.log('FROZEN_WORKFLOW_METADATA_VIEW ' + JSON.stringify({ workflow, revision,
      live_workflow_sha256: liveWorkflowSha, view_workflow_sha256: sha(workflowBytes),
      pinned_input_views: pinnedInputViews,
      scope: 'workflow metadata and explicitly pinned historical inputs are substituted in an isolated view; current checkout and evidence remain unchanged' }));
    // The enclosing workflow's existing job budget governs the validator.
    const result = spawnSync(command[0], command.slice(1), { cwd: view, stdio: 'inherit',
      env: { ...process.env, PYTHONDONTWRITEBYTECODE: '1' } });
    if (result.error) throw result.error;
    if (result.status === null) throw new Error(`validator terminated: ${result.signal}`);
    return result.status;
  } finally {
    rmSync(view, { recursive: true, force: true });
  }
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const args = process.argv.slice(2);
  const options = { roots: [] };
  try {
    while (args.length && args[0] !== '--') {
      const flag = args.shift();
      if (!['--workflow', '--revision', '--manifest', '--input-view-manifest', '--root'].includes(flag) || !args.length) {
        throw new Error(`unknown or incomplete argument: ${flag}`);
      }
      const value = args.shift();
      if (flag === '--root') options.roots.push(value);
      else if (flag === '--input-view-manifest') options.inputViewManifest = value;
      else options[flag.slice(2)] = value;
    }
    if (args.shift() !== '--') throw new Error('expected -- before the validator command');
    options.command = args;
    process.exitCode = withFrozenWorkflow(options);
  } catch (error) {
    console.error(`frozen workflow custody failed: ${error.message}`);
    process.exitCode = 1;
  }
}
