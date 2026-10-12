import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { withFrozenWorkflow } from './frozen-workflow-check.mjs';

const workflow = '.github/workflows/study.yml';
const study = 'research/study';
const manifest = study + '/FROZEN.json';
const historicalWorkflow = 'name: frozen historical workflow\n';
const liveWorkflow = 'name: current scoped workflow\n';
const liveSource = 'CURRENT_INTERFACE_V2\n';
const sha = value => createHash('sha256').update(value).digest('hex');

async function fixture(run) {
  const cwd = mkdtempSync(path.join(os.tmpdir(), 'crypto-frozen-fixture-'));
  const write = (file, text) => {
    const target = path.join(cwd, file);
    mkdirSync(path.dirname(target), { recursive: true });
    writeFileSync(target, text);
  };
  const git = (...args) => execFileSync('git', args, { cwd, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }).trim();
  try {
    git('init', '-q');
    write(workflow, historicalWorkflow);
    write('src/current.rs', 'HISTORICAL_INTERFACE_V1\n');
    write(study + '/input.json', '{"evidence":"unchanged"}\n');
    write(study + '/check.mjs', [
      "import assert from 'node:assert/strict';",
      "import { readFileSync, writeFileSync } from 'node:fs';",
      'assert.equal(readFileSync(' + JSON.stringify(workflow) + ", 'utf8'), " + JSON.stringify(historicalWorkflow) + ');',
      "assert.equal(readFileSync('src/current.rs', 'utf8'), " + JSON.stringify(liveSource) + ');',
      'writeFileSync(process.argv[2], process.cwd());',
    ].join('\n') + '\n');
    const files = Object.fromEntries(['input.json', 'check.mjs'].map(file => [file, sha(readFileSync(path.join(cwd, study, file)))]));
    files[workflow] = sha(historicalWorkflow);
    write(manifest, JSON.stringify({ files }) + '\n');
    git('add', '.');
    git('-c', 'user.name=CI fixture', '-c', 'user.email=ci-fixture@example.invalid', 'commit', '-qm', 'original freeze');
    const revision = git('rev-parse', 'HEAD');
    write(workflow, liveWorkflow);
    write('src/current.rs', liveSource);
    git('add', '.');
    git('-c', 'user.name=CI fixture', '-c', 'user.email=ci-fixture@example.invalid', 'commit', '-qm', 'current interface and CI metadata');
    const current = git('rev-parse', 'HEAD');
    const marker = path.join(cwd, 'view-location.txt');
    const options = { cwd, workflow, revision, manifest, roots: [study, 'src'],
      command: [process.execPath, study + '/check.mjs', marker] };
    await run({ cwd, write, options, current, marker });
  } finally {
    rmSync(cwd, { recursive: true, force: true });
  }
}

test('historical metadata is isolated while current source and live checkout are retained', async () => {
  await fixture(async ({ cwd, options, marker }) => {
    assert.equal(await withFrozenWorkflow(options), 0);
    const view = readFileSync(marker, 'utf8');
    assert.notEqual(view, cwd);
    assert.equal(existsSync(view), false, 'temporary source view must be cleaned up');
    assert.equal(readFileSync(path.join(cwd, workflow), 'utf8'), liveWorkflow);
    assert.equal(readFileSync(path.join(cwd, 'src/current.rs'), 'utf8'), liveSource);
  });
});

test('a validator failure remains a failed status and cleans its source view', async () => {
  await fixture(async ({ options, marker }) => {
    const code = "require('node:fs').writeFileSync(process.argv[1], process.cwd()); process.exit(17)";
    assert.equal(await withFrozenWorkflow({ ...options, command: [process.execPath, '-e', code, marker] }), 17);
    assert.equal(existsSync(readFileSync(marker, 'utf8')), false);
  });
});

test('workflow CLI arguments reach the isolated validator', async () => {
  await fixture(async ({ cwd, options, marker }) => {
    const helper = fileURLToPath(new URL('./frozen-workflow-check.mjs', import.meta.url));
    execFileSync(process.execPath, [helper, '--workflow', workflow, '--revision', options.revision,
      '--manifest', manifest, '--root', study, '--root', 'src', '--', ...options.command],
    { cwd, stdio: ['ignore', 'pipe', 'pipe'] });
    assert.equal(existsSync(readFileSync(marker, 'utf8')), false);
    assert.equal(readFileSync(path.join(cwd, workflow), 'utf8'), liveWorkflow);
  });
});

test('a shallow checkout fetches its exact missing pin and keeps current source', async () => {
  await fixture(async ({ cwd, options, marker }) => {
    const parent = mkdtempSync(path.join(os.tmpdir(), 'crypto-frozen-shallow-'));
    const clone = path.join(parent, 'repo');
    try {
      execFileSync('git', ['clone', '-q', '--depth=1', '--no-tags', pathToFileURL(cwd).href, clone],
        { stdio: ['ignore', 'pipe', 'pipe'] });
      assert.throws(() => execFileSync('git', ['cat-file', '-e', options.revision + '^{commit}'],
        { cwd: clone, stdio: ['ignore', 'pipe', 'pipe'] }));
      assert.equal(await withFrozenWorkflow({ ...options, cwd: clone }), 0);
      assert.equal(existsSync(readFileSync(marker, 'utf8')), false);
      assert.equal(readFileSync(path.join(clone, workflow), 'utf8'), liveWorkflow);
      assert.equal(readFileSync(path.join(clone, 'src/current.rs'), 'utf8'), liveSource);
    } finally {
      rmSync(parent, { recursive: true, force: true });
    }
  });
});

test('coordinated manifest and frozen-source replacement fails pinned manifest custody', async () => {
  await fixture(async ({ cwd, write, options, marker }) => {
    const replacement = '{"evidence":"replacement"}\n';
    write(study + '/input.json', replacement);
    const freeze = JSON.parse(readFileSync(path.join(cwd, manifest), 'utf8'));
    freeze.files['input.json'] = sha(replacement);
    write(manifest, JSON.stringify(freeze) + '\n');
    await assert.rejects(async () => withFrozenWorkflow(options));
    assert.equal(existsSync(marker), false);
  });
});

test('changed frozen evidence is rejected without changing the live workflow', async () => {
  await fixture(async ({ cwd, write, options, marker }) => {
    write(study + '/input.json', '{"evidence":"tampered"}\n');
    await assert.rejects(async () => withFrozenWorkflow(options));
    assert.equal(existsSync(marker), false);
    assert.equal(readFileSync(path.join(cwd, workflow), 'utf8'), liveWorkflow);
  });
});

test('deleted frozen evidence fails closed', async () => {
  await fixture(async ({ cwd, options, marker }) => {
    rmSync(path.join(cwd, study, 'input.json'));
    await assert.rejects(async () => withFrozenWorkflow(options));
    assert.equal(existsSync(marker), false);
  });
});

test('deleted manifest fails closed', async () => {
  await fixture(async ({ cwd, options, marker }) => {
    rmSync(path.join(cwd, manifest));
    await assert.rejects(async () => withFrozenWorkflow(options));
    assert.equal(existsSync(marker), false);
  });
});

test('a revision with mismatched workflow bytes cannot satisfy the freeze', async () => {
  await fixture(async ({ options, current, marker }) => {
    await assert.rejects(async () => withFrozenWorkflow({ ...options, revision: current }));
    assert.equal(existsSync(marker), false);
  });
});

test('paths outside the source view are rejected', async () => {
  await fixture(async ({ options, marker }) => {
    await assert.rejects(async () => withFrozenWorkflow({ ...options, roots: ['../'] }));
    assert.equal(existsSync(marker), false);
  });
});
