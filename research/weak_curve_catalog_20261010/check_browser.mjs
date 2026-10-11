#!/usr/bin/env node
// Dependency-free browser checks for presentation only. Never dispatch research.
import { spawn } from 'node:child_process';
import { existsSync, mkdtempSync, mkdirSync, readFileSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { resolve, dirname } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import assert from 'node:assert/strict';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const chrome = process.env.WCF_CATALOG_CHROME || [
  '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  '/usr/bin/google-chrome', '/usr/bin/chromium', '/usr/bin/chromium-browser'
].find(existsSync);
assert.ok(chrome, 'Set WCF_CATALOG_CHROME to an installed Chrome/Chromium executable');
const output = process.env.WCF_CATALOG_SCREENSHOTS;
if (output) mkdirSync(output, { recursive: true });
const profile = mkdtempSync(resolve(tmpdir(), 'weak-catalog-browser-'));
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
  const url=pathToFileURL(resolve(root,'docs/curves/weak-families/index.html')).href;
  await call('Page.navigate',{url});
  for (let i=0;i<100 && !(await evaluate("document.readyState === 'complete' && !!document.getElementById('stats')?.textContent"));i++) await sleep(50);
  await frame();
  assert.equal(await evaluate("document.querySelectorAll('.card').length"),21);
  assert.ok(await evaluate("document.getElementById('stats').textContent.includes('814 exact models')"));
  await screenshot('desktop-families');
  await evaluate("$('search').value='norm-one';$('search').dispatchEvent(new Event('input'));true");
  assert.ok(await evaluate("document.querySelectorAll('.card').length>=2"));
  await evaluate("$('search').value='';$('filter').value='cover_family';$('filter').dispatchEvent(new Event('change'));true");
  assert.equal(await evaluate("document.querySelectorAll('.card').length"),9);
  await evaluate("$('curvesTab').click();true");
  assert.ok(await evaluate("$('count').textContent.startsWith('814 ')"));
  await evaluate("$('next').click();true");
  assert.ok(await evaluate("$('page').textContent.startsWith('Page 2 ')"));
  await evaluate("$('filter').value='independent_model_law_large_source';$('filter').dispatchEvent(new Event('change'));true");
  assert.ok(await evaluate("$('count').textContent.startsWith('512 ')"));
  await evaluate("document.querySelector('.model').click();true");
  assert.ok(await evaluate("$('curveDetail').textContent.includes('unspecified / Unmeasured')"));
  assert.ok(await evaluate("$('curveDetail').textContent.includes('Trace:')"));
  await screenshot('desktop-model-detail');
  await evaluate("$('classesTab').click();true");
  assert.equal(await evaluate("document.querySelectorAll('#classRows tr').length"),294);
  await evaluate("$('filter').value='zero';$('filter').dispatchEvent(new Event('change'));true");
  assert.equal(await evaluate("document.querySelectorAll('#classRows tr').length"),46);
  await screenshot('desktop-class-zeros');
  await evaluate("$('filter').value='both';$('filter').dispatchEvent(new Event('change'));true");
  assert.equal(await evaluate("document.querySelectorAll('#classRows tr').length"),88);
  receipts.push('desktop: 21 family cards, search, 9 covering entries, 814 models, 512 source filter, pagination and exact model detail; 294 classes with 46 zeros and 88 overlaps');
  await evaluate("$('familiesTab').click();window.scrollTo(0,0);true");
  for (const width of [390,320]) {
    await call('Emulation.setDeviceMetricsOverride',{width,height:1000,deviceScaleFactor:1,mobile:true});
    await frame();
    assert.ok(await evaluate("document.documentElement.scrollWidth<=innerWidth+1"),'page overflow at '+width);
    await screenshot('mobile-'+width);
  }
  receipts.push('390px and 320px mobile: cards readable with no page overflow');
  await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
  await call('Emulation.setEmulatedMedia',{features:[{name:'prefers-color-scheme',value:'dark'}]});
  await screenshot('desktop-dark');
  receipts.push('dark theme rendered');
  await call('Emulation.setScriptExecutionDisabled',{value:true});
  await call('Page.navigate',{url});
  for(let i=0;i<100&&!(await evaluate("document.readyState==='complete'"));i++)await sleep(50);
  assert.ok(await evaluate("document.querySelector('noscript').textContent.includes('Filtering requires JavaScript')"));
  assert.ok(await evaluate("!!document.querySelector('a[href=\"models.json\"]')"));
  receipts.push('JavaScript disabled: native source links and scope remain available');
  assert.equal(errors.length,0,'Browser script errors: '+JSON.stringify(errors));
  assert.equal(network.filter(url=>/^https?:/.test(url)).length,0,'Catalogue made external requests');
  const summary={status:'PASS_WEAK_CURVE_CATALOG_BROWSER',chrome,checks:receipts,
    script_errors:errors,external_requests:[],viewport_widths:[1280,390,320],
    page_sha256:(await import('node:crypto')).createHash('sha256').update(readFileSync(resolve(root,'docs/curves/weak-families/index.html'))).digest('hex')};
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
