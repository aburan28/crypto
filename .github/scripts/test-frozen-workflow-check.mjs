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
const historicalLock = 'released nested Cargo.lock\n';
const currentLock = 'current root Cargo.lock\n';
const lockPath = 'research/archive/Cargo.lock';
const originalFreeze = 'research/rotation/FROZEN.json';
const inputViewManifest = '.github/scripts/n53-historical-cargo-lock-view-v1.json';
const sha = value => createHash('sha256').update(value).digest('hex');

async function fixture(run, { withPinnedInput = false } = {}) {
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
    write('Cargo.lock', historicalLock);
    write(lockPath, historicalLock);
    const frozenSourceBytes = JSON.stringify({ source_sha256: { cargo_lock: sha(historicalLock) } }) + '\n';
    write(originalFreeze, frozenSourceBytes);
    write(study + '/input.json', '{"evidence":"unchanged"}\n');
    write(study + '/check.mjs', [
      "import assert from 'node:assert/strict';",
      "import { readFileSync, writeFileSync } from 'node:fs';",
      'assert.equal(readFileSync(' + JSON.stringify(workflow) + ", 'utf8'), " + JSON.stringify(historicalWorkflow) + ');',
      "assert.equal(readFileSync('src/current.rs', 'utf8'), " + JSON.stringify(liveSource) + ');',
      'assert.equal(readFileSync(' + JSON.stringify(lockPath) + ", 'utf8'), "
        + JSON.stringify(withPinnedInput ? historicalLock : currentLock) + ');',
      "assert.equal(readFileSync('Cargo.lock', 'utf8'), " + JSON.stringify(currentLock) + ');',
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
    write('Cargo.lock', currentLock);
    write(lockPath, currentLock);
    write(inputViewManifest, JSON.stringify({
      schema: 'frozen-input-view/v1',
      inputs: [{
        target: lockPath,
        revision,
        source_path: lockPath,
        sha256: sha(historicalLock),
        original_manifest: originalFreeze,
        original_manifest_sha256: sha(frozenSourceBytes),
        original_source_key: 'cargo_lock',
        current_equal: [[lockPath, 'Cargo.lock']],
      }],
    }) + '\n');
    git('add', '.');
    git('-c', 'user.name=CI fixture', '-c', 'user.email=ci-fixture@example.invalid', 'commit', '-qm', 'current interface and CI metadata');
    const current = git('rev-parse', 'HEAD');
    const marker = path.join(cwd, 'view-location.txt');
    const options = { cwd, workflow, revision, manifest, roots: [study, 'src',
      'research/archive', 'research/rotation', 'Cargo.lock'],
      command: [process.execPath, study + '/check.mjs', marker] };
    if (withPinnedInput) options.inputViewManifest = inputViewManifest;
    await run({ cwd, write, options, current, marker, withPinnedInput });
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
      '--manifest', manifest, '--root', study, '--root', 'src', '--root', 'research/archive',
      '--root', 'research/rotation', '--root', 'Cargo.lock', '--', ...options.command],
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

test('released input overlay is isolated and records current lock provenance', async () => {
  await fixture(async ({ cwd, options, marker }) => {
    assert.equal(await withFrozenWorkflow(options), 0);
    const view = readFileSync(marker, 'utf8');
    assert.notEqual(view, cwd);
    assert.equal(existsSync(view), false, 'temporary historical input view must be cleaned up');
    assert.equal(readFileSync(path.join(cwd, lockPath), 'utf8'), currentLock);
    assert.equal(readFileSync(path.join(cwd, 'Cargo.lock'), 'utf8'), currentLock);
  }, { withPinnedInput: true });
});

test('historical Git input must match the original freeze digest', async () => {
  await fixture(async ({ cwd, write, options, marker, current }) => {
    const config = JSON.parse(readFileSync(path.join(cwd, options.inputViewManifest), 'utf8'));
    config.inputs[0].revision = current;
    write(options.inputViewManifest, JSON.stringify(config) + '\n');
    await assert.rejects(async () => withFrozenWorkflow(options), /pinned Git source hash mismatch/);
    assert.equal(existsSync(marker), false);
    assert.equal(readFileSync(path.join(cwd, lockPath), 'utf8'), currentLock);
  }, { withPinnedInput: true });
});

test('current nested/root lock drift fails before the historical overlay', async () => {
  await fixture(async ({ cwd, write, options, marker }) => {
    write('Cargo.lock', 'different current root lock\n');
    await assert.rejects(async () => withFrozenWorkflow(options), /current source inputs differ/);
    assert.equal(existsSync(marker), false);
    assert.equal(readFileSync(path.join(cwd, lockPath), 'utf8'), currentLock);
  }, { withPinnedInput: true });
});

test('current lock comparison cannot be reduced to a self-comparison', async () => {
  await fixture(async ({ cwd, write, options, marker }) => {
    const config = JSON.parse(readFileSync(path.join(cwd, options.inputViewManifest), 'utf8'));
    config.inputs[0].current_equal = [[lockPath, lockPath]];
    write(options.inputViewManifest, JSON.stringify(config) + '\n');
    await assert.rejects(async () => withFrozenWorkflow(options), /invalid current-source equality pair/);
    assert.equal(existsSync(marker), false);
  }, { withPinnedInput: true });
});

test('pinned input digest cannot diverge from the original frozen manifest', async () => {
  await fixture(async ({ cwd, write, options, marker }) => {
    const config = JSON.parse(readFileSync(path.join(cwd, options.inputViewManifest), 'utf8'));
    config.inputs[0].sha256 = '0'.repeat(64);
    write(options.inputViewManifest, JSON.stringify(config) + '\n');
    await assert.rejects(async () => withFrozenWorkflow(options), /does not match original frozen source/);
    assert.equal(existsSync(marker), false);
  }, { withPinnedInput: true });
});

test('a changed original source manifest rejects the pinned input view', async () => {
  await fixture(async ({ cwd, write, options, marker }) => {
    write(originalFreeze, '{"source_sha256":{"cargo_lock":"tampered"}}\n');
    await assert.rejects(async () => withFrozenWorkflow(options), /original freeze manifest changed/);
    assert.equal(existsSync(marker), false);
  }, { withPinnedInput: true });
});
