#!/usr/bin/env node
// Local, NUL-delimited git diffs avoid the API's 300-file path-filter limit.
// Missing comparisons select every check; they must never silently skip work.
import { execFileSync } from 'node:child_process';
import { appendFileSync, existsSync, readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

export const scopes = ['rust_root', 'rust_cd', 'rust_fmt', 'release_rust',
  'release_native', 'release_cuda', 'ecbench_replay', 'actions', 'c_cpp',
  'javascript', 'python', 'rust'];

const maintained = /^(src|tests|examples|benches)\//;
const frozenExamples = /^experiments\/koblitz-s3-pair-query-20261007-v4\/(baseline|candidate)\.rs$/;
const rustConfig = /^(Cargo\.(toml|lock)|build\.rs|rust-toolchain(\.toml)?|rustfmt\.toml|clippy\.toml|\.cargo\/.*)$/;
// The release packages this independent crate and its path dependency. These
// paths do not belong to the root Cargo package, so classify their own sources,
// manifests, build script and embedded inputs explicitly.
const isogenyCliRust = /^tools\/isogeny-cli\//;
const isogenyAlgosInputs = /^isogeny_algos\/(?:Cargo\.toml|build\.rs|src\/.*\.rs)$/;
const isogenyEmbeddedInputs = new Set([
  'docs/curves/registry.json',
  'docs/curves/standards/coverage.json',
  'research/p192_p224_hundred_large_prime_isogenies_20261010/PROTOCOL.md',
  'research/p192_p224_large_degree_isogenies_20261009/PROTOCOL.md',
  'research/p192_p224_large_degree_isogenies_20261009/verification_poly.rs',
  'research/p192_p224_large_degree_isogenies_20261009/verification_kernel.rs',
]);
const prose = /\.(md|rst|txt|pdf|svg|png|jpg|jpeg|gif)$/i;
const control = /^\.github\/scripts\/.*ci-scope.*$/;
// These trees are copied into every Rust release archive. Documentation is
// normally compilation-neutral, but package payload changes (and attributes
// that govern their byte identity) must exercise cross-platform staging and
// the installed-tree smoke path before merge.
const rustReleasePayload = /^docs\/(?:curves|ecdlp-audit)\//;
const rustReleaseRootPayload = /^(?:README\.md|LICENSE(?:\.md)?)$/;

export function classify(files, embedded = new Set(), cargoSources = new Set()) {
  const result = Object.fromEntries(scopes.map(key => [key, false]));
  for (const file of files) {
    if (control.test(file)) return Object.fromEntries(scopes.map(key => [key, true]));
    const releaseAttributes = file === '.gitattributes';
    const root = /\.rs$/.test(file) && maintained.test(file) || frozenExamples.test(file) || cargoSources.has(file) || rustConfig.test(file) || embedded.has(file);
    result.rust_root ||= root || file === 'scripts/clippy-maintained.sh';
    result.rust_cd ||= /^cd_attack\/(src|tests|examples|benches)\//.test(file) && !prose.test(file)
      || /^cd_attack\/(Cargo\.(toml|lock)|build\.rs|rust-toolchain(\.toml)?|\.cargo\/.*)$/.test(file)
      || /^rust-toolchain(\.toml)?$|^\.cargo\//.test(file);
    result.rust_fmt ||= /\.rs$/.test(file) && (maintained.test(file) || /^cd_attack\//.test(file) || file === 'build.rs');
    result.release_rust ||= /^src\//.test(file) && !prose.test(file) || rustConfig.test(file)
      || file === 'tests/fixtures/ecbench_manifest_small.json'
      || rustReleasePayload.test(file) || rustReleaseRootPayload.test(file)
      || releaseAttributes || embedded.has(file)
      || isogenyCliRust.test(file) || isogenyAlgosInputs.test(file)
      || isogenyEmbeddedInputs.has(file);
    // Both Makefile targets compile src/main.cu against include/ and the
    // committed generated headers; research benchmark sources are separate.
    const client = file === 'ecc2k130/src/main.cu' || file === 'ecc2k130/Makefile'
      || /^ecc2k130\/(include|generated)\/.*\.(h|hpp|cuh|inc)$/.test(file);
    result.release_native ||= client || releaseAttributes || file === 'ecc2k130/README.md'
      || /^gpu\/(ecc|ecc2k|semaev|macaulay|btcpuzzle)\//.test(file) && !prose.test(file);
    result.release_cuda ||= client || releaseAttributes;
    result.ecbench_replay ||= root
      || /^docs\/(ecbench|bounds)\//.test(file) && !prose.test(file)
      || /^research\/(ecbench_[^/]+|ghs_curve_comparison_20261009)\//.test(file) && (!prose.test(file) || file.includes('/sessions/'))
      || /^research\/notes\/ecc2k130\/(cold_compact_orbit_ic_20261003|full_rank_compact_orbit_ic_20261003)\//.test(file) && (!prose.test(file) || file.includes('/sessions/'));
    result.actions ||= /^\.github\/(workflows|actions)\//.test(file);
    result.c_cpp ||= /\.(c|cc|cpp|cxx|h|hh|hpp|hxx|cu|cuh)$/.test(file)
      || /(^|\/)(Makefile|CMakeLists\.txt)$|\.cmake$/.test(file);
    result.javascript ||= /\.(js|mjs|cjs|jsx|ts|tsx)$/.test(file)
      || /(^|\/)(package(-lock)?\.json|pnpm-lock\.yaml|yarn\.lock|bun\.lockb?)$/.test(file);
    result.python ||= /\.py$/.test(file)
      || /(^|\/)(pyproject\.toml|requirements[^/]*\.txt|uv\.lock|poetry\.lock|setup\.(cfg|py)|tox\.ini|\.python-version)$/.test(file);
    result.rust ||= /\.rs$/.test(file) || /(^|\/)Cargo\.(toml|lock)$/.test(file) || rustConfig.test(file);
    if (/^\.github\/workflows\/rust-lint\.ya?ml$/.test(file)) {
      result.rust_root = result.rust_cd = result.rust_fmt = true;
    }
    if (/^\.github\/workflows\/rust-all-targets\.ya?ml$/.test(file)) result.rust_root = true;
    if (/^\.github\/workflows\/release\.ya?ml$/.test(file)) {
      result.release_rust = result.release_native = result.release_cuda = true;
    }
    if (/^\.github\/workflows\/ecbench\.ya?ml$/.test(file)) result.ecbench_replay = true;
    if (/^\.github\/workflows\/codeql\.ya?ml$/.test(file)) {
      result.actions = result.c_cpp = result.javascript = result.python = result.rust = true;
    }
  }
  return result;
}

function git(args, cwd) {
  return execFileSync('git', args, { cwd, encoding: 'utf8', maxBuffer: 16 * 1024 * 1024,
    env: { ...process.env, GIT_TERMINAL_PROMPT: '0' }, stdio: ['ignore', 'pipe', 'pipe'] });
}

function trackedFiles(cwd) {
  return git(['ls-files', '-z'], cwd).split('\0').filter(Boolean);
}

function sourceText(cwd, file) {
  return existsSync(path.join(cwd, file)) ? readFileSync(path.join(cwd, file), 'utf8') : git(['show', `HEAD:${file}`], cwd);
}

// Cargo paths are literal TOML strings. Read the manifest from HEAD when a
// scope-only sparse checkout hasn't materialized it, without loading archives.
export function cargoSourceFiles(cwd, tracked = trackedFiles(cwd)) {
  if (!tracked.includes('Cargo.toml')) return new Set();
  const manifest = existsSync(path.join(cwd, 'Cargo.toml'))
    ? readFileSync(path.join(cwd, 'Cargo.toml'), 'utf8')
    : git(['show', 'HEAD:Cargo.toml'], cwd);
  const targets = new Set();
  for (const block of manifest.split(/^\s*(?:\[\[(?:example|bin|bench|test)\]\]|\[lib\])\s*$/m).slice(1)) {
    const table = block.split(/^\s*\[/m)[0];
    const line = table.match(/^\s*path\s*=\s*(.*?)\s*$/m)?.[1];
    if (line === undefined) continue;
    const literal = line.match(/^("(?:[^"\\]|\\.)*"|'[^']*')\s*(?:#.*)?$/);
    if (!literal) throw new Error('unsupported Cargo target path syntax');
    const value = literal[1][0] === '"' ? JSON.parse(literal[1]) : literal[1].slice(1, -1);
    const target = path.posix.normalize(value);
    if (target.startsWith('../') || path.posix.isAbsolute(target)) throw new Error('Cargo target escapes the checkout');
    targets.add(target);
  }
  const sources = new Set(targets);
  // Out-of-tree targets may declare sibling/nested modules, such as the
  // frozen n53 examples' reference.rs. Fetch source files, never run archives.
  for (const target of targets) {
    const dir = path.posix.dirname(target);
    for (const file of tracked) {
      if (/\.rs$/.test(file) && (dir === '.' ? path.posix.dirname(file) === '.' : file.startsWith(`${dir}/`))) sources.add(file);
    }
  }
  return sources;
}

function materializeSources(cwd, sources) {
  const missing = [...sources].filter(file => !existsSync(path.join(cwd, file)));
  if (missing.length) git(['sparse-checkout', 'add', ...missing.map(file => `/${file}`)], cwd);
}

export function sourceFiles(cwd) {
  const tracked = trackedFiles(cwd);
  const trackedSet = new Set(tracked);
  const sources = new Set(tracked.filter(file => /\.rs$/.test(file)
    && (maintained.test(file) || frozenExamples.test(file) || file === 'build.rs')));
  for (const source of cargoSourceFiles(cwd, tracked)) sources.add(source);
  const queue = [...sources];
  for (let index = 0; index < queue.length; index++) {
    const file = queue[index];
    const text = sourceText(cwd, file);
    const dir = path.posix.dirname(file);
    const add = dependency => {
      const normalized = path.posix.normalize(dependency);
      if (trackedSet.has(normalized) && !sources.has(normalized)) {
        sources.add(normalized);
        queue.push(normalized);
      }
    };
    const requireSource = candidates => {
      const found = candidates.map(candidate => path.posix.normalize(candidate)).filter(candidate => trackedSet.has(candidate));
      if (!found.length) throw new Error(`missing Rust source dependency in ${file}: ${candidates[0]}`);
      for (const dependency of found) add(dependency);
    };
    // Follow exact external modules recursively. The alternate module base
    // also covers path attributes within an inline module in foo.rs.
    for (const match of text.matchAll(/#\[\s*path\s*=\s*"([^"\n]+)"\s*\]/g)) {
      requireSource([path.posix.join(dir, match[1]), path.posix.join(dir, path.posix.basename(file, '.rs'), match[1])]);
    }
    for (const match of text.matchAll(/include!\(\s*"([^"\n]+)"/g)) requireSource([path.posix.join(dir, match[1])]);
    for (const match of text.matchAll(/include!\(\s*concat!\(\s*env!\("CARGO_MANIFEST_DIR"\),\s*"([^"\n]+)"/g)) requireSource([match[1].replace(/^\//, '')]);
    // Ordinary modules reached outside the maintained source directories
    // need their own siblings/submodules, too. Only tracked source candidates
    // are included; historical run outputs never enter this graph.
    const overrides = new Set([...text.matchAll(/#\[\s*path\s*=\s*"[^"\n]+"\s*\](?:\s|\/\/[^\n]*\n|#\[[^\]]*\])*(?:pub(?:\([^)]*\))?\s+)?mod\s+([A-Za-z_][A-Za-z_0-9]*)\s*;/g)].map(match => match[1]));
    for (const match of text.matchAll(/^\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+([A-Za-z_][A-Za-z_0-9]*)\s*;/gm)) {
      if (overrides.has(match[1])) continue;
      const candidates = [dir, path.posix.join(dir, path.posix.basename(file, '.rs'))]
        .flatMap(base => [path.posix.join(base, `${match[1]}.rs`), path.posix.join(base, match[1], 'mod.rs')]);
      requireSource(candidates);
    }
  }
  return sources;
}

// Compile-time include inputs are dependencies even when they are JSON,
// SQL or historical research fixtures. Keep them derived from source rather
// than a hand-maintained exclusion list that can miss the next include_str!.
export function embeddedInputs(cwd, sources = sourceFiles(cwd)) {
  const inputs = new Set();
  for (const file of sources) {
    const text = sourceText(cwd, file);
    for (const match of text.matchAll(/include_(?:str|bytes)!\(\s*"([^"\n]+)"/g)) {
      inputs.add(path.posix.normalize(path.posix.join(path.posix.dirname(file), match[1])));
    }
    for (const match of text.matchAll(/include_(?:str|bytes)!\(\s*concat!\(\s*env!\("CARGO_MANIFEST_DIR"\),\s*"([^"\n]+)"/g)) {
      inputs.add(match[1].replace(/^\//, ''));
    }
  }
  return inputs;
}

export function materializeInputs(cwd = process.cwd()) {
  const sources = sourceFiles(cwd);
  materializeSources(cwd, sources);
  const inputs = [...embeddedInputs(cwd, sources)];
  if (inputs.some(input => input.startsWith('../') || path.posix.isAbsolute(input))) {
    throw new Error('compile-time input escapes the checkout');
  }
  if (inputs.length) git(['sparse-checkout', 'add', ...inputs.map(input => `/${input}`)], cwd);
  console.log(`materialized ${inputs.length} compile-time inputs`);
}

export function select({ event, base, head = 'HEAD', cwd = process.cwd() }) {
  const all = Object.fromEntries(scopes.map(key => [key, true]));
  if (!['pull_request', 'push', 'merge_group'].includes(event)) return all;
  if (!/^[0-9a-f]{40,64}$/i.test(base ?? '') || /^0+$/.test(base)) return all;
  try {
    try { git(['cat-file', '-e', `${base}^{commit}`], cwd); }
    catch { git(['fetch', '--no-tags', '--depth=1', 'origin', base], cwd); }
    // For PRs, the workflow passes the tested merge commit's first parent as
    // base. A direct diff works with its shallow, depth-two checkout.
    const files = git(['diff', '--name-only', '-z', '--no-renames', base, head, '--'], cwd)
      .split('\0').filter(Boolean);
    const sources = sourceFiles(cwd);
    materializeSources(cwd, sources);
    return classify(files, embeddedInputs(cwd, sources), sources);
  } catch (error) {
    console.error(`::warning::CI scope comparison unavailable; running all checks (${error.message.split('\n')[0]})`);
    return all;
  }
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  if (process.argv.includes('--materialize')) {
    materializeInputs();
    process.exit(0);
  }
  const result = select({ event: process.env.CI_EVENT, base: process.env.CI_BASE });
  const text = Object.entries(result).map(([key, value]) => `${key}=${value}\n`).join('');
  if (process.env.GITHUB_OUTPUT) appendFileSync(process.env.GITHUB_OUTPUT, text);
  process.stdout.write(text);
}
