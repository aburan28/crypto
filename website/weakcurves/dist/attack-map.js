(() => {
  const $ = selector => document.querySelector(selector);
  const grid = $('#attack-grid');
  const search = $('#attack-search');
  const family = $('#attack-family');
  const state = $('#attack-state');
  const more = $('#attack-more');
  const criterionFamilies = {
    'Inputs and subgroups': ['smooth_selected_order','proper_subfield_points','pairing_subgroup_validation'],
    'Auxiliary information': ['static_dh_oracle'],
    'Pairing target fields': ['mov_weil_pairing','frey_ruck_tate_pairing','supersingular_pairing','ddh_pairing','minimal_pairing_field','pairing_group_roles','pairing_target_field_dlp'],
    'Index calculus and covers': ['ghs_binary','ghs_odd_covering','cubic_norm_one','quadratic_torus','extension_index_calculus','cover_and_decomposition','cubic_jacobian_isogeny'],
    'Endomorphism acceleration': ['glv_endomorphism','gls_frobenius','isogeny_transfer'],
    'Intrinsic DLP': ['anomalous_prime'],
    'Isogeny protocols': ['sidh_auxiliary_torsion']
  };
  const familyOf = id => Object.entries(criterionFamilies).find(([,ids]) => ids.includes(id))?.[0] || 'Other';
  let items = [], shown = 16;

  function add(parent, tag, value, className) {
    const element = document.createElement(tag);
    element.textContent = value;
    if (className) element.className = className;
    parent.append(element);
    return element;
  }

  function card(item) {
    const article = document.createElement('article');
    article.className = 'attack-card';
    const heading = document.createElement('div');
    heading.className = 'attack-card-head';
    add(heading,'span',item.family,'attack-family');
    add(heading,'span',item.state === 'indexed' ? 'Criterion indexed' : 'Research queue',`attack-state ${item.state}`);
    article.append(heading);
    add(article,'h3',item.name);
    add(article,'p',item.scope,'attack-scope');
    add(article,'p',item.requires,'attack-requires');
    const detail = document.createElement('details');
    add(detail,'summary','Evidence gate and source ↓');
    add(detail,'p',`Test: ${item.test}`);
    if (item.construction) add(detail,'p',`Construction: ${item.construction}`);
    const source = add(detail,'a',`${item.source} ↗`);
    source.href=item.url;source.target='_blank';source.rel='noopener noreferrer';
    if (item.state === 'indexed') {
      const jump = add(detail,'a','Open indexed criterion ↗','attack-jump');
      jump.href=`#criterion-${item.id}`;
    }
    article.append(detail);
    return article;
  }

  function render() {
    const q = search.value.trim().toLowerCase();
    const filtered = items.filter(item => (family.value === 'all' || item.family === family.value)
      && (state.value === 'all' || item.state === state.value)
      && (!q || [item.name,item.id,item.family,item.scope,item.requires,item.test,item.source]
        .some(value => String(value).toLowerCase().includes(q))));
    grid.replaceChildren(...filtered.slice(0,shown).map(card));
    $('#attack-count').textContent=`${filtered.length} matching checks · ${items.filter(i=>i.state==='indexed').length} in current criteria · ${items.filter(i=>i.state==='queue').length} queued for model-level assessment`;
    more.hidden=shown>=filtered.length;
    if(!filtered.length) add(grid,'p','No attack check matches these filters.');
  }

  Promise.all([
    fetch('criteria.json').then(response => {if(!response.ok) throw Error('criteria unavailable');return response.json();}),
    fetch('attack_map.json').then(response => {if(!response.ok) throw Error('attack map unavailable');return response.json();})
  ]).then(([criteria,expansion]) => {
    items = [
      ...criteria.map(c => ({id:c.id,name:c.name,family:familyOf(c.id),scope:c.domain,
        requires:c.recognition,test:c.witness,construction:c.construction,source:c.source,url:c.url,state:'indexed'})),
      ...expansion.items.map(item => ({...item,state:'queue'}))
    ];
    $('#attack-total').textContent=items.length;
    for (const name of [...new Set(items.map(item=>item.family))].sort()) {
      const option=document.createElement('option');option.value=name;option.textContent=name;family.append(option);
    }
    render();
  }).catch(error => {$('#attack-count').textContent=`Could not load attack map: ${error.message}`;});
  [search,family,state].forEach(control => control.addEventListener('input',() => {shown=16;render();}));
  more.addEventListener('click',() => {shown+=16;render();});
})();
