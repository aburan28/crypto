import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { existsSync, mkdtempSync, mkdirSync, writeFileSync, rmSync, renameSync } from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import { classify, embeddedInputs, materializeInputs, scopes, select } from './ci-scope.mjs';

test('prose does not run compilers, profiling or CodeQL', () => {
  for (const file of ['README.md', 'AGENTS.md', 'docs/ecbench/README.md',
    'examples/README.md', 'research/ecbench_example/REPORT.md']) {
    assert.ok(Object.values(classify([file])).every(value => !value), file);
  }
});

test('crates and languages select their own checks', () => {
  const cd = classify(['cd_attack/src/lib.rs']);
  assert.equal(cd.rust_cd, true);
  assert.equal(cd.rust_root, false);
  assert.equal(cd.release_cuda, false);
  assert.equal(classify(['research/archive/native/src/lib.rs']).rust_root, false);
  const source = classify(['src/hash/sha256.rs']);
  assert.equal(source.rust_root, true);
  assert.equal(source.rust_cd, false);
  assert.equal(source.release_rust, true);
  assert.equal(source.release_native, false);
  assert.equal(source.rust, true);
  assert.equal(source.python, false);
  assert.equal(classify(['tools/export.py']).python, true);
  assert.equal(classify(['tools/dashboard.mjs']).javascript, true);
  assert.equal(classify(['ecc2k130/Makefile']).release_cuda, true);
  assert.equal(classify(['ecc2k130/aws/README.md']).release_native, false);
});

test('dependencies, workflow edits and frozen evidence remain applicable', () => {
  assert.equal(classify(['Cargo.lock']).rust_root, true);
  assert.equal(classify(['build.rs']).rust_root, true);
  assert.equal(classify(['build.rs']).rust_fmt, true);
  assert.equal(classify(['.cargo/config.toml']).rust_root, true);
  assert.equal(classify(['.cargo/config.toml']).rust_cd, true);
  assert.equal(classify(['.github/workflows/rust-all-targets.yml']).rust_root, true);
  assert.equal(classify(['.github/workflows/other.yml']).actions, true);
  assert.equal(classify(['.github/workflows/other.yml']).rust_root, false);
  assert.equal(classify(['docs/bounds/records/record.json']).ecbench_replay, true);
  assert.equal(classify(['research/ecbench_example/sessions/one/records.jsonl']).ecbench_replay, true);
  assert.equal(classify(['research/ecbench_example/sessions/one/README.md']).ecbench_replay, true);
  assert.equal(classify(['research/fixture.json'], new Set(['research/fixture.json'])).rust_root, true);
});

test('real git diffs include renames, deletes, unusual names and over 300 files', () => {
  const cwd = mkdtempSync(path.join(os.tmpdir(), 'crypto-ci-scope-'));
  const git = (...args) => execFileSync('git', args, { cwd, encoding: 'utf8' }).trim();
  const write = (file, text) => { mkdirSync(path.dirname(path.join(cwd, file)), { recursive: true }); writeFileSync(path.join(cwd, file), text); };
  try {
    git('init', '-q'); git('config', 'user.email', 'ci@example.invalid'); git('config', 'user.name', 'CI scope test');
    write('src/lib.rs', 'const DATA: &str = include_str!("../research/fixture.json");\n');
    write('src/concat.rs', 'const DATA: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/docs/input.json"));\n');
    write('research/fixture.json', '{}\n'); write('tools/old.py', 'print(1)\n');
    git('add', '.'); git('commit', '-qm', 'base');
    const base = git('rev-parse', 'HEAD');
    assert.deepEqual([...embeddedInputs(cwd)].sort(), ['docs/input.json', 'research/fixture.json']);
    renameSync(path.join(cwd, 'tools/old.py'), path.join(cwd, 'tools/new.py'));
    write('docs/line\nbreak.md', 'doc\n');
    for (let n = 0; n < 350; n++) write(`docs/file-${n}.md`, 'doc\n');
    write('research/fixture.json', '{"value":1}\n');
    git('add', '.'); git('commit', '-qm', 'change');
    const result = select({ event: 'pull_request', base, cwd });
    assert.equal(result.python, true); assert.equal(result.rust_root, true);
    assert.equal(result.release_cuda, false);
    rmSync(path.join(cwd, 'src/lib.rs'));
    git('add', '.'); git('commit', '-qm', 'delete');
    assert.equal(select({ event: 'merge_group', base, cwd }).rust_root, true);
    for (const args of [{ event: 'schedule' }, { event: 'push', base: '0'.repeat(40) },
      { event: 'pull_request', base: 'f'.repeat(40) }, { event: 'pull_request', base, head: 'missing-ref' }]) {
      assert.deepEqual(select({ ...args, cwd }), Object.fromEntries(scopes.map(key => [key, true])));
    }
  } finally { rmSync(cwd, { recursive: true, force: true }); }
});

test('sparse compilation restores only embedded inputs, preserving archived files off disk', () => {
  const cwd = mkdtempSync(path.join(os.tmpdir(), 'crypto-ci-inputs-'));
  const git = (...args) => execFileSync('git', args, { cwd, encoding: 'utf8' }).trim();
  const write = (file, text) => { mkdirSync(path.dirname(path.join(cwd, file)), { recursive: true }); writeFileSync(path.join(cwd, file), text); };
  try {
    git('init', '-q'); git('config', 'user.email', 'ci@example.invalid'); git('config', 'user.name', 'CI scope test');
    write('src/lib.rs', 'const DATA: &str = include_str!("../research/fixture.json");\n');
    write('build.rs', 'const DATA: &str = include_str!("research/build.json");\n');
    write('research/build.json', '{}\n');
    write('research/fixture.json', '{}\n'); write('research/archive/huge.json', '{}\n');
    git('add', '.'); git('commit', '-qm', 'base');
    git('sparse-checkout', 'set', '--no-cone', '/src/', '/build.rs');
    assert.equal(existsSync(path.join(cwd, 'research/fixture.json')), false);
    materializeInputs(cwd);
    assert.equal(existsSync(path.join(cwd, 'research/fixture.json')), true);
    assert.equal(existsSync(path.join(cwd, 'research/build.json')), true);
    assert.equal(existsSync(path.join(cwd, 'research/archive/huge.json')), false);
  } finally { rmSync(cwd, { recursive: true, force: true }); }
});
