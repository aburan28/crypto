// Offline behavioral checks for the generated performance dashboard.
// Run: node scripts/site/test_performance_ui.js
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const path = require('node:path');
const html = fs.readFileSync(path.join(__dirname, '../../docs/performance-gains.html'), 'utf8');
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
function fixture(href, clipboard) {
  const nodes = {}, listeners = {};
  const window = {location: new URL(href), addEventListener(name, fn) { listeners[name] = fn; }};
  window.history = {replaceState(_state, _title, next) { window.location = new URL(next); }};
  const document = {addEventListener() {}, getElementById(id) {
    return nodes[id] ??= {value: id === 'local-metric' ? 'specialization_s' : id === 'full-metric' ? 'wall_s' : '', checked: id === 'intervals', hidden: true, innerHTML: '', textContent: '', style: {},
      events: {}, setAttribute() {}, focus() {}, select() {},
      addEventListener(name, fn) { (this.events[name] ??= []).push(fn); }};
  }};
  const ctx = vm.createContext({window, document, URL, navigator: {clipboard}});
  vm.runInContext(script, ctx);
  return {nodes, ctx, window, listeners,
    async trigger(id, event) { for (const fn of nodes[id].events[event] ?? []) await fn(); }};
}
(async () => {
  const f = fixture('https://example.test/crypto/scoreboard/performance-gains.html?local=setup_s&full=cpu_s&legacy=1&pairs=1&intervals=0#confirmation');
  assert.equal(f.nodes['local-metric'].value, 'setup_s');
  assert.equal(f.nodes['full-metric'].value, 'cpu_s');
  assert.equal(f.nodes.legacy.checked, true);
  assert.equal(f.nodes.intervals.checked, false);
  assert.equal(f.nodes['raw-pairs'].checked, true);
  for (const metric of ['specialization_s', 'setup_s']) for (const legacy of [true, false]) for (const intervals of [true, false]) {
    f.nodes['local-metric'].value = metric; f.nodes.legacy.checked = legacy; f.nodes.intervals.checked = intervals;
    await f.trigger('local-metric', 'change');
    assert.doesNotMatch(f.nodes['local-chart'].innerHTML, /NaN|undefined|Infinity/);
  }
  for (const metric of ['wall_s', 'cpu_s']) for (const raw of [true, false]) {
    f.nodes['full-metric'].value = metric; f.nodes['raw-pairs'].checked = raw;
    await f.trigger('full-metric', 'change');
    assert.doesNotMatch(f.nodes['full-chart'].innerHTML, /NaN|undefined|Infinity/);
  }
  const roundtrip = fixture(f.window.location.href);
  assert.equal(roundtrip.nodes['full-metric'].value, f.nodes['full-metric'].value);
  assert.equal(roundtrip.nodes['raw-pairs'].checked, f.nodes['raw-pairs'].checked);
  assert.equal(f.window.location.hash, '#confirmation');
  f.window.location = new URL('https://example.test/?local=invalid&full=invalid&intervals=invalid');
  f.listeners.popstate();
  assert.equal(f.nodes['local-metric'].value, 'specialization_s');
  assert.equal(f.nodes['full-metric'].value, 'wall_s');
  assert.equal(f.nodes.intervals.checked, true);
  await f.trigger('share-view', 'click');
  assert.equal(f.nodes['share-fallback'].hidden, false);
  assert.match(f.nodes['share-status'].textContent, /Copy the link below/);
  let copied;
  const success = fixture('https://example.test/performance.html', {writeText: async url => { copied = url; }});
  await success.trigger('share-view', 'click');
  assert.match(copied, /local=specialization_s/);
  assert.match(success.nodes['share-status'].textContent, /Link copied/);
  const denied = fixture('https://example.test/performance.html', {writeText: async () => { throw Error('denied'); }});
  await denied.trigger('share-view', 'click');
  assert.equal(denied.nodes['share-fallback'].hidden, false);
  const local = fixture('file:///tmp/performance-gains.html');
  await local.trigger('share-view', 'click');
  assert.match(local.nodes['share-status'].textContent, /Local preview/);
  assert.equal(local.nodes['share-fallback'].hidden, false);
  const ids = [...html.matchAll(/\bid="([^"]+)"/g)].map(x => x[1]);
  assert.equal(new Set(ids).size, ids.length);
  for (const anchor of html.matchAll(/href="#([^"]+)"/g)) assert(ids.includes(anchor[1]), anchor[1]);
  assert.doesNotMatch(html, /__LOCAL_TABLE__|__FULL_TABLE__|__FROZEN_DATA__/);
  console.log('PASS: 12 chart states, URL round-trip, invalid parameters, history, clipboard success/denial, offline fallback, section anchors.');
})().catch(error => { console.error(error); process.exitCode = 1; });
