import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { existsSync, mkdtempSync, mkdirSync, readFileSync, writeFileSync, rmSync, renameSync } from 'node:fs';
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
  for (const file of ['ecc2k130/src/main.cu', 'ecc2k130/include/packedengine.cuh', 'ecc2k130/generated/eccF131.h']) {
    assert.equal(classify([file]).release_native, true, file);
    assert.equal(classify([file]).release_cuda, true, file);
  }
  assert.equal(classify(['ecc2k130/benchmarks/old/probe.cu']).release_cuda, false);
  assert.equal(classify(['ecc2k130/aws/README.md']).release_native, false);
});

test('dependencies, workflow edits and frozen evidence remain applicable', () => {
  assert.equal(classify(['Cargo.lock']).rust_root, true);
  assert.equal(classify(['build.rs']).rust_root, true);
  assert.equal(classify(['build.rs']).rust_fmt, true);
  assert.equal(classify(['scripts/clippy-maintained.sh']).rust_root, true);
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

test('isogeny release inputs select the Rust package lane', () => {
  const inputs = [
    'tools/isogeny-cli/Cargo.toml',
    'tools/isogeny-cli/Cargo.lock',
    'tools/isogeny-cli/build.rs',
    'tools/isogeny-cli/src/main.rs',
    'tools/isogeny-cli/README.md',
    'isogeny_algos/Cargo.toml',
    'isogeny_algos/src/lib.rs',
    'isogeny_algos/src/field.rs',
    'docs/curves/registry.json',
    'docs/curves/standards/coverage.json',
    'research/p192_p224_hundred_large_prime_isogenies_20261010/PROTOCOL.md',
    'research/p192_p224_large_degree_isogenies_20261009/PROTOCOL.md',
    'research/p192_p224_large_degree_isogenies_20261009/verification_poly.rs',
    'research/p192_p224_large_degree_isogenies_20261009/verification_kernel.rs',
  ];
  for (const file of inputs) {
    const result = classify([file]);
    assert.equal(result.release_rust, true, file);
    assert.equal(result.release_native, false, file);
    assert.equal(result.release_cuda, false, file);
  }
  // Unrelated isogeny research reports do not invoke the release builder.
  assert.equal(classify(['research/p192_p224_large_degree_isogenies_20261008/PROTOCOL.md']).release_rust, false);
});

test('volcano-cli, a path dependency of the root crate, selects the root lanes', () => {
  for (const file of ['tools/volcano-cli/Cargo.toml', 'tools/volcano-cli/Cargo.lock', 'tools/volcano-cli/src/main.rs']) {
    const result = classify([file]);
    assert.equal(result.rust_root, true, file);
    assert.equal(result.release_rust, true, file);
    assert.equal(result.ecbench_replay, true, file);
    assert.equal(result.release_native, false, file);
    assert.equal(result.release_cuda, false, file);
  }
  assert.equal(classify(['isogeny_algos/src/lib.rs']).rust_root, true);
  assert.equal(classify(['tools/volcano-cli/README.md']).rust_root, false);
  assert.equal(classify(['tools/volcano-cli/README.md']).release_rust, true);
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

test('new Cargo targets and source modules survive sparse checkouts without loading archives', () => {
  const cwd = mkdtempSync(path.join(os.tmpdir(), 'crypto-ci-targets-'));
  const git = (...args) => execFileSync('git', args, { cwd, encoding: 'utf8' }).trim();
  const write = (file, text) => { mkdirSync(path.dirname(path.join(cwd, file)), { recursive: true }); writeFileSync(path.join(cwd, file), text); };
  try {
    git('init', '-q'); git('config', 'user.email', 'ci@example.invalid'); git('config', 'user.name', 'CI scope test');
    write('Cargo.toml', '[package]\nname="fixture"\nversion="0.1.0"\n[[example]]\nname="fresh"\npath = "experiments/fresh/native.rs" # target\n');
    write('src/lib.rs', '// maintained\n');
    write('experiments/fresh/native.rs', 'mod reference;\n');
    write('experiments/fresh/reference.rs', 'const DATA: &str = include_str!("input.json");\n');
    write('experiments/fresh/input.json', '{}\n');
    write('experiments/fresh/archive/raw.jsonl', '{}\n');
    git('add', '.'); git('commit', '-qm', 'base');
    const base = git('rev-parse', 'HEAD');
    write('experiments/fresh/reference.rs', 'const DATA: &str = include_str!("input.json");\n// revised\n');
    git('add', '.'); git('commit', '-qm', 'module change');
    // Scope-only jobs can read Cargo.toml from Git without its checkout.
    git('sparse-checkout', 'set', '--no-cone', '/src/');
    assert.equal(select({ event: 'pull_request', base, cwd }).rust_root, true);
    assert.equal(existsSync(path.join(cwd, 'experiments/fresh/native.rs')), true);
    assert.equal(existsSync(path.join(cwd, 'experiments/fresh/reference.rs')), true);
    assert.equal(existsSync(path.join(cwd, 'experiments/fresh/input.json')), false);
    materializeInputs(cwd);
    assert.equal(existsSync(path.join(cwd, 'experiments/fresh/input.json')), true);
    assert.equal(existsSync(path.join(cwd, 'experiments/fresh/archive/raw.jsonl')), false);
  } finally { rmSync(cwd, { recursive: true, force: true }); }
});

test('external Rust path modules and includes are followed recursively', () => {
  const cwd = mkdtempSync(path.join(os.tmpdir(), 'crypto-ci-modules-'));
  const git = (...args) => execFileSync('git', args, { cwd, encoding: 'utf8' }).trim();
  const write = (file, text) => { mkdirSync(path.dirname(path.join(cwd, file)), { recursive: true }); writeFileSync(path.join(cwd, file), text); };
  try {
    git('init', '-q'); git('config', 'user.email', 'ci@example.invalid'); git('config', 'user.name', 'CI scope test');
    write('examples/check.rs', '#[path="../research/study/reference.rs"]\nmod reference;\n');
    write('research/study/reference.rs', 'mod arithmetic;\ninclude!("render.inc");\n');
    write('research/study/arithmetic.rs', 'const DATA: &str = include_str!("fixture.json");\n');
    write('research/study/render.inc', '#[path="display.rs"]\nmod display;\n');
    write('research/study/display.rs', '// display\n');
    write('research/study/fixture.json', '{}\n');
    write('research/study/archive/raw.jsonl', '{}\n');
    git('add', '.'); git('commit', '-qm', 'base');
    const base = git('rev-parse', 'HEAD');
    write('research/study/display.rs', '// revised display\n');
    git('add', '.'); git('commit', '-qm', 'external module change');
    git('sparse-checkout', 'set', '--no-cone', '/examples/');
    assert.equal(select({ event: 'pull_request', base, cwd }).rust_root, true);
    for (const file of ['reference.rs', 'arithmetic.rs', 'render.inc', 'display.rs']) {
      assert.equal(existsSync(path.join(cwd, 'research/study', file)), true, file);
    }
    materializeInputs(cwd);
    assert.equal(existsSync(path.join(cwd, 'research/study/fixture.json')), true);
    assert.equal(existsSync(path.join(cwd, 'research/study/archive/raw.jsonl')), false);
    // Deleting an external path module, include! file or ordinary submodule
    // must fail closed even when its declaring source remains unchanged.
    for (const file of ['reference.rs', 'render.inc', 'arithmetic.rs']) {
      const relative = `research/study/${file}`;
      const saved = readFileSync(path.join(cwd, relative), 'utf8');
      const deletionBase = git('rev-parse', 'HEAD');
      rmSync(path.join(cwd, relative));
      git('add', '.'); git('commit', '-qm', `delete ${file}`);
      assert.deepEqual(select({ event: 'pull_request', base: deletionBase, cwd }), Object.fromEntries(scopes.map(key => [key, true])));
      write(relative, saved); git('add', '.'); git('commit', '-qm', `restore ${file}`);
    }
  } finally { rmSync(cwd, { recursive: true, force: true }); }
});
