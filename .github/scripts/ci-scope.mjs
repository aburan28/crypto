#!/usr/bin/env node
// Local, NUL-delimited git diffs avoid the API's 300-file path-filter limit.
// Missing comparisons select every check; they must never silently skip work.
import { execFileSync } from 'node:child_process';
import { appendFileSync, readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

export const scopes = ['rust_root', 'rust_cd', 'rust_fmt', 'release_rust',
  'release_native', 'release_cuda', 'ecbench_replay', 'actions', 'c_cpp',
  'javascript', 'python', 'rust'];

const maintained = /^(src|tests|examples|benches)\//;
const frozenExamples = /^experiments\/koblitz-s3-pair-query-20261007-v4\/(baseline|candidate)\.rs$/;
const rustConfig = /^(Cargo\.(toml|lock)|build\.rs|rust-toolchain(\.toml)?|rustfmt\.toml|clippy\.toml|\.cargo\/.*)$/;
const prose = /\.(md|rst|txt|pdf|svg|png|jpg|jpeg|gif)$/i;
const control = /^\.github\/scripts\/.*ci-scope.*$/;

export function classify(files, embedded = new Set()) {
  const result = Object.fromEntries(scopes.map(key => [key, false]));
  for (const file of files) {
    if (control.test(file)) return Object.fromEntries(scopes.map(key => [key, true]));
    const root = /\.rs$/.test(file) && maintained.test(file) || frozenExamples.test(file) || rustConfig.test(file) || embedded.has(file);
    result.rust_root ||= root;
    result.rust_cd ||= /^cd_attack\/(src|tests|examples|benches)\//.test(file) && !prose.test(file)
      || /^cd_attack\/(Cargo\.(toml|lock)|build\.rs|rust-toolchain(\.toml)?|\.cargo\/.*)$/.test(file)
      || /^rust-toolchain(\.toml)?$|^\.cargo\//.test(file);
    result.rust_fmt ||= /\.rs$/.test(file) && (maintained.test(file) || /^cd_attack\//.test(file) || file === 'build.rs');
    result.release_rust ||= /^src\//.test(file) && !prose.test(file) || rustConfig.test(file)
      || file === 'tests/fixtures/ecbench_manifest_small.json'
      || /^docs\/curves\//.test(file) || embedded.has(file);
    result.release_native ||= /^ecc2k130\/[^/]+\.(c|cc|cpp|h|hpp|cu|cuh)$/.test(file) || file === 'ecc2k130/Makefile'
      || /^gpu\/(ecc|ecc2k|semaev|macaulay|btcpuzzle)\//.test(file) && !prose.test(file);
    result.release_cuda ||= /^ecc2k130\/[^/]+\.(c|cc|cpp|h|hpp|cu|cuh)$/.test(file) || file === 'ecc2k130/Makefile';
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

// Compile-time include inputs are dependencies even when they are JSON,
// SQL or historical research fixtures. Keep them derived from source rather
// than a hand-maintained exclusion list that can miss the next include_str!.
export function embeddedInputs(cwd) {
  const inputs = new Set();
  const sources = git(['ls-files', '-z'], cwd).split('\0').filter(file =>
    /\.rs$/.test(file) && (maintained.test(file) || frozenExamples.test(file) || file === 'build.rs'));
  for (const file of sources) {
    const text = readFileSync(path.join(cwd, file), 'utf8');
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
  const inputs = [...embeddedInputs(cwd)];
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
    // HEAD is the tested merge result for PRs; comparing it directly to the
    // event's base also works in shallow checkouts without a merge-base walk.
    const files = git(['diff', '--name-only', '-z', '--no-renames', base, head, '--'], cwd)
      .split('\0').filter(Boolean);
    return classify(files, embeddedInputs(cwd));
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
