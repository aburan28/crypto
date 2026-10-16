// Presentation regression only: preserve existing scientific panels across regeneration.
import {readFileSync} from 'node:fs';
import {execFileSync} from 'node:child_process';
import assert from 'node:assert/strict';
const pagePath='docs/index-calculus-scoreboard.html', dataPath='docs/ic/dashboard-overview-data.json';
const before=readFileSync(pagePath,'utf8'), data=readFileSync(dataPath,'utf8');
const publication=execFileSync('git',['show','14e56ddb0902161eaa289425c5131a58fefd5243:'+pagePath],{maxBuffer:4*1024*1024}).toString();
function panel(text,id){const match=text.match(new RegExp('<section\\b[^>]*id="'+id+'"[^>]*>[\\s\\S]*?</section>'));assert.ok(match,'Missing '+id);return match[0];}
for(const id of ['n83-matched-20261014','n83-basis-setup-20261015'])assert.equal(panel(before,id),panel(publication,id),'Published evidence panel changed');
execFileSync('python3',['scripts/render_n83_basis_diagnostics.py'],{stdio:'pipe'});
execFileSync('node',['tools/render_ic_dashboard.mjs'],{stdio:'pipe'});
execFileSync('node',['tools/render_ic_dashboard.mjs','--check'],{stdio:'pipe'});
const after=readFileSync(pagePath,'utf8');assert.equal(after,before,'Presenters do not compose idempotently');assert.equal(readFileSync(dataPath,'utf8'),data,'Generated overview data drifted');
for(const id of ['n83-matched-20261014','n83-basis-setup-20261015'])assert.equal(after.split('id="'+id+'"').length-1,1,'Duplicate or missing N83 panel');
console.log('PASS: both N83 evidence panels unchanged, present once, regeneration pipeline idempotent; source pins and historical ledger retained');
