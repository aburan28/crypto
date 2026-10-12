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

export function withFrozenWorkflow({ cwd = process.cwd(), workflow, revision, manifest, roots = [], command }) {
  workflow = relative(workflow);
  manifest = relative(manifest);
  roots = roots.map(relative);
  if (!/^[0-9a-f]{40,64}$/.test(revision ?? '')) throw new Error('expected an immutable Git revision');
  if (!Array.isArray(command) || !command.length || command.some(arg => typeof arg !== 'string')) {
    throw new Error('expected a command and argument array');
  }
  try {
    execFileSync('git', ['cat-file', '-e', `${revision}^{commit}`], { cwd, stdio: ['ignore', 'pipe', 'pipe'] });
  } catch {
    execFileSync('git', ['fetch', '--filter=blob:none', '--no-tags', '--depth=1', 'origin', revision],
      { cwd, stdio: ['ignore', 'pipe', 'pipe'] });
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
    // Original preflights retain their Git ancestry assertions, against the
    // current checkout's unchanged object database; no historical code checkout.
    const gitPath = path.join(cwd, '.git');
    if (!existsSync(gitPath)) throw new Error('current Git metadata is missing');
    symlinkSync(gitPath, path.join(view, '.git'));
    console.log('FROZEN_WORKFLOW_METADATA_VIEW ' + JSON.stringify({ workflow, revision,
      live_workflow_sha256: liveWorkflowSha, view_workflow_sha256: sha(workflowBytes),
      scope: 'only workflow scheduling metadata substituted; current source/evidence copied unchanged; original replay workflow digest refers to this view' }));
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
      if (!['--workflow', '--revision', '--manifest', '--root'].includes(flag) || !args.length) {
        throw new Error(`unknown or incomplete argument: ${flag}`);
      }
      const value = args.shift();
      if (flag === '--root') options.roots.push(value);
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
