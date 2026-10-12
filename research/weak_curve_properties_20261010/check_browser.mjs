#!/usr/bin/env node
// Dependency-free browser checks for presentation only. Never dispatch research.
import { spawn } from 'node:child_process';
import { existsSync, mkdtempSync, mkdirSync, readFileSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { resolve, dirname } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import assert from 'node:assert/strict';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const chrome = process.env.WCF_PROPERTIES_CHROME || [
  '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  '/usr/bin/google-chrome', '/usr/bin/chromium', '/usr/bin/chromium-browser'
].find(existsSync);
assert.ok(chrome, 'Set WCF_PROPERTIES_CHROME to an installed Chrome/Chromium executable');
const output = process.env.WCF_PROPERTIES_SCREENSHOTS;
if (output) mkdirSync(output, { recursive: true });
const profile = mkdtempSync(resolve(tmpdir(), 'weak-properties-browser-'));
const child = spawn(chrome, ['--headless=new', '--remote-debugging-port=0',
  '--remote-debugging-address=127.0.0.1', '--no-first-run', '--no-default-browser-check',
  '--disable-background-networking', '--disable-component-update',
  '--disable-extensions', '--user-data-dir=' + profile, 'about:blank'],
  { stdio: ['ignore','ignore','pipe'] });
let chromeLog = '', socket;
child.stderr.on('data', bytes => { chromeLog += bytes; });
const sleep = ms => new Promise(r => setTimeout(r,ms));
const deadline = Date.now() + 20000;
const receipts = [], errors = [], network = [];
try {
  while (!/DevTools listening on (ws:\/\/[^\s]+)/.test(chromeLog)) {
    assert.ok(child.exitCode === null, 'Chrome failed to start: ' + chromeLog);
    assert.ok(Date.now() < deadline, 'Chrome launch deadline exceeded: ' + chromeLog);
    await sleep(50);
  }
  const debuggerUrl = new URL(chromeLog.match(/DevTools listening on (ws:\/\/[^\s]+)/)[1]);
  const target = await (await fetch(`http://${debuggerUrl.host}/json/new?about:blank`,
    { method: 'PUT', signal: AbortSignal.timeout(5000) })).json();
  socket = new WebSocket(target.webSocketDebuggerUrl);
  await new Promise((resolve,reject) => { socket.onopen=resolve; socket.onerror=reject; });
  let nextId = 0;
  const pending = new Map();
  socket.onmessage = event => {
    const data = JSON.parse(event.data);
    if (data.id) {
      const call = pending.get(data.id);
      if (!call) return;
      pending.delete(data.id); clearTimeout(call.timer);
      if (data.error) call.reject(new Error(JSON.stringify(data.error))); else call.resolve(data.result);
    }
    if (data.method === 'Runtime.exceptionThrown') errors.push(data.params.exceptionDetails);
    if (data.method === 'Network.requestWillBeSent') network.push(data.params.request.url);
  };
  function call(method,params={}) {
    return new Promise((resolve,reject) => {
      const id=++nextId;
      const timer=setTimeout(()=>{pending.delete(id);reject(new Error('CDP deadline: '+method));},10000);
      pending.set(id,{resolve,reject,timer}); socket.send(JSON.stringify({id,method,params}));
    });
  }
  async function evaluate(expression) {
    const result=await call('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});
    assert.ok(!result.exceptionDetails, JSON.stringify(result.exceptionDetails));
    return result.result.value;
  }
  async function frame() { await evaluate('new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))'); }
  async function screenshot(name) {
    if (!output) return;
    await frame();
    const {data}=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
    writeFileSync(resolve(output,name+'.png'),Buffer.from(data,'base64'));
  }
  await call('Page.enable'); await call('Runtime.enable'); await call('Network.enable');
  await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
  await call('Emulation.setEmulatedMedia',{features:[{name:'prefers-color-scheme',value:'light'},{name:'prefers-reduced-motion',value:'reduce'}]});
  const html=resolve(root,'docs/curves/weak-families/properties.html');
  const url=pathToFileURL(html).href;
  await call('Page.navigate',{url});
  for(let i=0;i<100&&!(await evaluate("document.readyState==='complete' && !!$('stats')?.textContent"));i++)await sleep(50);
  await frame();
  async function filter(id,value) { await evaluate(`$('${id}').value=${JSON.stringify(value)};$('${id}').dispatchEvent(new Event('${id==='query'?'input':'change'}'));true`); }
  async function count(n) { assert.ok(await evaluate(`$('count').textContent.startsWith('${n} model records;')`),'expected '+n+' records'); }
  assert.ok(await evaluate("$('stats').textContent.startsWith('814 exact model records')"));
  await count(814);
  assert.equal(await evaluate("document.querySelectorAll('#rows tr').length"),50);
  await screenshot('desktop-properties');
  await evaluate("$('next').click();true");
  assert.ok(await evaluate("$('page').textContent.startsWith('Page 2 /')"));
  await filter('property','exactlevel'); await count(360);
  await filter('property','floor'); await count(102);
  await filter('property','full4'); await count(117);
  await filter('cohort','independent_large'); await count(65);
  await screenshot('desktop-large-full4');
  await filter('property','exactlevel'); await count(219);
  await filter('property',''); await count(512);
  await filter('depth','9'); await count(1);
  await filter('depth',''); await filter('cohort','');
  await filter('property','jproper'); await count(4);
  await filter('property','literalproper'); await count(1);
  await filter('property','isoproper'); await count(3);
  await filter('property','jblocked'); await count(1);
  await evaluate("document.querySelector('.model').click();true");
  assert.ok(await evaluate("$('detail').textContent.includes('icv1-fp3k6-t202-39e02f24')"));
  assert.ok(await evaluate("$('detail').textContent.includes('j 3, displayed coefficients 6, isomorphic model 6')"));
  assert.ok(await evaluate("$('detail').textContent.includes('permitted 2-depth from 0 to 0')"));
  await screenshot('desktop-twist-descent');
  await filter('property','smooth'); await count(10);
  await filter('property','ss'); await count(0);
  assert.ok(await evaluate("$('previous').disabled && $('next').disabled"));
  await filter('property','floor'); await filter('query','higher'); await count(0);
  await filter('query',''); await filter('cohort','higher_degree_control'); await count(25);
  await evaluate("document.querySelector('.model').click();true");
  assert.ok(await evaluate("/permitted 2-depth from ([1-9]) to \\1/.test($('detail').textContent)"));
  receipts.push('814 records; 50-row pagination; 360 certified levels; 102 floors; 117 full-4 models; 512 large sources with 65 full-4 and 219 certified levels; one depth-9 source');
  receipts.push('field-definition filters: 4 proper j, 1 proper displayed model, 3 proper isomorphic models, one 3/6/6 twist obstruction; 10 13-smooth orders; zero supersingular sources; positive floor bounds rendered correctly');
  await filter('cohort',''); await filter('property',''); await filter('splitting','inert'); await count(99);
  await filter('splitting',''); await filter('depth','ge2'); await count(539);
  await filter('depth','');
  await evaluate("$('detail').classList.add('hide');window.scrollTo(0,0);true");
  for(const width of [390,320]) {
    await call('Emulation.setDeviceMetricsOverride',{width,height:1000,deviceScaleFactor:1,mobile:true});
    await frame();
    assert.ok(await evaluate("document.documentElement.scrollWidth<=innerWidth+1"),'page overflow at '+width);
    assert.ok(await evaluate("document.querySelector('.scroll').scrollWidth>document.querySelector('.scroll').clientWidth"),'table should scroll internally');
    await screenshot('mobile-'+width);
  }
  receipts.push('390px and 320px mobile: no page overflow; wide property table scrolls within its container');
  await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
  await call('Emulation.setEmulatedMedia',{features:[{name:'prefers-color-scheme',value:'dark'}]});
  await screenshot('desktop-dark'); receipts.push('dark theme rendered');
  await call('Emulation.setScriptExecutionDisabled',{value:true});
  await call('Page.navigate',{url});
  for(let i=0;i<100&&!(await evaluate("document.readyState==='complete'"));i++)await sleep(50);
  assert.ok(await evaluate("document.querySelector('noscript').textContent.includes('Filtering requires JavaScript')"));
  assert.ok(await evaluate("!!document.querySelector('a[href=\"properties.json\"]')"));
  receipts.push('JavaScript disabled: scope and canonical source links remain available');
  assert.equal(errors.length,0,'Browser script errors: '+JSON.stringify(errors));
  assert.equal(network.filter(url=>/^https?:/.test(url)).length,0,'Property browser made external requests');
  const summary={status:'PASS_WEAK_CURVE_PROPERTIES_BROWSER',chrome,checks:receipts,script_errors:errors,external_requests:[],viewport_widths:[1280,390,320],page_sha256:(await import('node:crypto')).createHash('sha256').update(readFileSync(html)).digest('hex')};
  if(output)writeFileSync(resolve(output,'browser-check.json'),JSON.stringify(summary,null,2)+'\n');
  console.log(JSON.stringify(summary,null,2));
} finally {
  socket?.close();
  child.kill('SIGTERM');
  for (let i=0;i<50 && child.exitCode===null && child.signalCode===null;i++) await sleep(100);
  if(child.exitCode===null && child.signalCode===null) child.kill('SIGKILL');
  // Only this invocation's fresh browser profile is removed.
  rmSync(profile,{recursive:true,force:true});
}
