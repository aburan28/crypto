const $ = (selector) => document.querySelector(selector);
const table = $('#atlas-table');
const tbody = table.querySelector('tbody');
const thead = table.querySelector('thead');
const search = $('#search');
const fieldFilter = $('#field-filter');
const evidenceFilter = $('#evidence-filter');
const more = $('#more');
const tabs = [...document.querySelectorAll('[data-tab]')];
let snapshot = null;
let tab = 'signals';
let shown = 50;

const kind = {
  anomalous_prime_field: {label:'Trace-one',grade:'Exact DLP criterion',badge:'verified',rank:0},
  smooth_selected_group: {label:'Smooth subgroup',grade:'Exact DLP criterion',badge:'verified',rank:0},
  exact_covering_model: {label:'GHS cover model',grade:'Covering geometry',badge:'geometry',rank:1},
  low_embedding_degree_candidate: {label:'MOV / FR candidate',grade:'Candidate',badge:'candidate',rank:2},
  unclassified_covering_model: {label:'Covering open',grade:'Open',badge:'open',rank:3},
  outside_covering_index: {label:'Outside index',grade:'Open',badge:'open',rank:3},
  no_indexed_covering_model: {label:'No indexed cover',grade:'Open',badge:'open',rank:4}
};

function textCell(row, value, className='') {
  const cell = document.createElement('td');
  cell.textContent = value;
  if (className) cell.className = className;
  row.append(cell);
  return cell;
}

function curveRecords() {
  return [...snapshot.prime_examples,...snapshot.models].sort((a,b) => {
    const priority = kind[a.signal].rank - kind[b.signal].rank;
    if (priority) return priority;
    return b.bits-a.bits || String(a.id).localeCompare(String(b.id));
  });
}

function labelFor(record) {
  if (record.a !== undefined) return `E / F${record.p}: y² = x³ + ${record.a}x + ${record.b}`;
  return record.slug;
}

function fieldFor(record) {
  return record.degree === 1 ? `F${record.p}` : `F(${record.p}^${record.degree})`;
}

function openRecord(record) {
  const dialog = $('#record-dialog');
  const spec = kind[record.signal];
  $('#dialog-title').textContent = labelFor(record);
  $('#dialog-grade').textContent = spec.grade;
  $('#dialog-grade').className = `badge ${spec.badge}`;
  const facts = [
    ['Identity',record.id],['Field',fieldFor(record)],['Field bits',String(record.bits)],
    ['Trace',String(record.trace)],['Curve order',String(record.order)],
    ['Criterion',spec.label]
  ];
  if (record.generator) facts.push(['Selected point',`(${record.generator.join(', ')})`]);
  if (record.generator_order) facts.push(['Point order',String(record.generator_order)]);
  if (record.group_exponent) facts.push(['Group exponent',String(record.group_exponent)]);
  if (record.group_structure) facts.push(['Group structure',record.group_structure.join(' × ')]);
  if (record.embedding_degree) facts.push(['Embedding degree',String(record.embedding_degree)]);
  if (record.cubic) facts.push(['Cubic model',record.cubic],['Quadratic model',record.quadratic]);
  if (record.cubic_screen) facts.push(['Cubic screen',record.cubic_screen]);
  if (record.quadratic_screen) facts.push(['Quadratic screen',record.quadratic_screen]);
  facts.push(['Subgroup route',record.subgroup_advantage]);
  if (record.sha256) facts.push(['Model SHA-256',record.sha256]);
  if (record.proof) facts.push(['Evidence',record.proof]);
  if (record.source) facts.push(['Source',record.source]);
  const dl = $('#dialog-facts');
  dl.replaceChildren();
  for (const [key,value] of facts) {
    const dt = document.createElement('dt');dt.textContent=key;
    const dd = document.createElement('dd');dd.textContent=value;
    dl.append(dt,dd);
  }
  $('#dialog-explain').textContent = record.signal === 'exact_covering_model'
    ? 'This exact model has a positive genus-three covering-model label. A selected subgroup map and comparative cost remain separate gates.'
    : record.signal === 'low_embedding_degree_candidate'
      ? 'The selected order and embedding degree are exact. The pairing image and destination-field cost remain to be certified.'
      : record.attack_ready
        ? 'The model and selected order satisfy a published discrete-log criterion. This snapshot reports no measured attack cost.'
        : 'This catalog record retains its exact identity and current classification status.';
  dialog.showModal();
}

function filteredCurves() {
  let rows = curveRecords();
  if (tab === 'signals') rows = rows.filter(r => kind[r.signal].rank <= 2);
  const q = search.value.trim().toLowerCase();
  if (q) rows = rows.filter(r => [r.id,r.slug,labelFor(r),fieldFor(r),kind[r.signal].label,r.trace].some(x => String(x||'').toLowerCase().includes(q)));
  const field = fieldFilter.value;
  if (field === 'prime') rows = rows.filter(r => r.degree === 1);
  if (field === 'extension') rows = rows.filter(r => r.degree > 1);
  if (field === 'large') rows = rows.filter(r => r.bits >= 192 && r.bits <= 252);
  const grade = evidenceFilter.value;
  if (grade === 'dlp') rows = rows.filter(r => kind[r.signal].rank === 0);
  if (grade === 'cover') rows = rows.filter(r => kind[r.signal].rank === 1);
  if (grade === 'candidate') rows = rows.filter(r => kind[r.signal].rank === 2);
  if (grade === 'open') rows = rows.filter(r => kind[r.signal].rank >= 3);
  return rows;
}

function renderCurves() {
  thead.replaceChildren();tbody.replaceChildren();
  const head = document.createElement('tr');
  ['#','Curve / identity','Field','Method','Evidence','Trace','Order'].forEach(v => {const th=document.createElement('th');th.textContent=v;head.append(th)});
  thead.append(head);
  const records = filteredCurves();
  records.slice(0,shown).forEach((r,i) => {
    const row = document.createElement('tr');
    textCell(row,String(i+1).padStart(2,'0'),'rank');
    const name = textCell(row,'','curve-name');
    const button = document.createElement('button');button.type='button';button.className='curve-button';
    button.textContent=labelFor(r);button.title='Open curve evidence';button.addEventListener('click',()=>openRecord(r));
    name.append(button);const small=document.createElement('small');small.textContent=r.id;name.append(small);
    textCell(row,fieldFor(r),'mono');
    textCell(row,kind[r.signal].label);
    const grade=textCell(row,'');const badge=document.createElement('span');badge.className=`badge ${kind[r.signal].badge}`;badge.textContent=kind[r.signal].grade;grade.append(badge);
    textCell(row,String(r.trace),'mono');textCell(row,String(r.order),'mono');
    tbody.append(row);
  });
  $('#result-count').textContent=`Showing ${Math.min(shown,records.length).toLocaleString()} of ${records.length.toLocaleString()} records`;
  more.hidden=shown>=records.length;
}

function renderCases() {
  thead.replaceChildren();tbody.replaceChildren();
  const head=document.createElement('tr');
  ['Case','Extension degree d','Cover rank n','Excess genus e','Cover genus','Coefficient'].forEach(v=>{const th=document.createElement('th');th.textContent=v;head.append(th)});
  thead.append(head);
  const q=search.value.trim().toLowerCase();
  const cases=snapshot.published_cases.filter(c=>!q||[`case ${c.case_no}`,c.degree,c.cover_genus,c.coefficient].some(v=>String(v).toLowerCase().includes(q)));
  cases.slice(0,shown).forEach(c=>{
    const tr=document.createElement('tr');
    const td=textCell(tr,'');const link=document.createElement('a');link.href=c.source_url;link.target='_blank';link.rel='noopener noreferrer';link.className='case-link';link.textContent=`Case ${String(c.case_no).padStart(2,'0')} ↗`;td.append(link);
    [c.degree,c.rank,c.excess_genus,c.cover_genus,c.coefficient==='eta'?'1 or nonsquare':'1 only'].forEach(v=>textCell(tr,String(v),'mono'));
    tbody.append(tr);
  });
  $('#result-count').textContent=`Showing ${Math.min(shown,cases.length)} of ${cases.length} published parameter cases`;
  more.hidden=shown>=cases.length;
}

function render() {
  if (!snapshot) return;
  $('#model-controls').style.display=tab==='families'?'none':'flex';
  if (tab==='families') renderCases();else renderCurves();
}

function setTab(next) {
  tab=next;shown=50;search.value='';fieldFilter.value='all';evidenceFilter.value='all';
  tabs.forEach(t=>{const active=t.dataset.tab===next;t.classList.toggle('active',active);t.setAttribute('aria-selected',String(active));});
  render();
}

function renderCriteria(criteria) {
  const grid=$('#criteria-grid');grid.replaceChildren();
  criteria.forEach((c,i)=>{
    const article=document.createElement('article');article.className='criterion';
    const top=document.createElement('div');top.className='criterion-top';
    const index=document.createElement('span');index.className='criterion-index';index.textContent=String(i+1).padStart(2,'0')+' / '+String(criteria.length).padStart(2,'0');
    const method=document.createElement('span');method.className='criterion-method';method.textContent=c.method;top.append(index,method);
    const h=document.createElement('h3');h.textContent=c.name;
    const domain=document.createElement('p');domain.className='domain';domain.textContent=c.domain;
    const recognition=document.createElement('p');recognition.textContent=c.recognition;
    const details=document.createElement('details');
    const summary=document.createElement('summary');summary.textContent='Witness & construction ↓';details.append(summary);
    const witness=document.createElement('p');witness.textContent='Witness: '+c.witness;
    const construction=document.createElement('p');construction.textContent='Construction: '+c.construction;
    const source=document.createElement('a');source.href=c.url;source.target='_blank';source.rel='noopener noreferrer';source.className='source-ref';source.textContent=c.source+' ↗';
    details.append(witness,construction,source);article.append(top,h,domain,recognition,details);grid.append(article);
  });
}

function initDialog() {
  const dialog=document.createElement('dialog');dialog.id='record-dialog';dialog.innerHTML='<div class="dialog-head"><div><p class="eyebrow">Curve evidence</p><h2 id="dialog-title"></h2></div><button type="button" class="close-dialog" aria-label="Close details">×</button></div><span id="dialog-grade"></span><p id="dialog-explain"></p><dl id="dialog-facts"></dl>';
  document.body.append(dialog);
  dialog.querySelector('.close-dialog').addEventListener('click',()=>dialog.close());
  dialog.addEventListener('click',event=>{if(event.target===dialog)dialog.close()});
}

Promise.all([fetch('snapshot.json').then(r=>{if(!r.ok)throw Error('snapshot unavailable');return r.json()}),fetch('criteria.json').then(r=>{if(!r.ok)throw Error('criteria unavailable');return r.json()})])
  .then(([data,criteria])=>{
    snapshot=data;initDialog();renderCriteria(criteria);
    $('#stat-models').textContent=data.counts.exact_models.toLocaleString();
    $('#stat-covers').textContent=data.counts.exact_covering_models.toLocaleString();
    $('#stat-cases').textContent=data.counts.published_cases.toLocaleString();
    $('#stat-large').textContent=data.counts.large_models.toLocaleString();
    $('#signal-count').textContent=data.counts.exact_covering_models+data.counts.prime_field_examples;
    $('#corpus-count').textContent=data.counts.exact_models+data.counts.prime_field_examples;
    $('#case-count').textContent=data.counts.published_cases;
    render();
  })
  .catch(error=>{$('#result-count').textContent=`Could not load atlas data: ${error.message}`;});

tabs.forEach(t=>t.addEventListener('click',()=>setTab(t.dataset.tab)));
[search,fieldFilter,evidenceFilter].forEach(el=>el.addEventListener('input',()=>{shown=50;render()}));
more.addEventListener('click',()=>{shown+=50;render()});
