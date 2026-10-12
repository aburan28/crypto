(() => {
  const $ = (selector) => document.querySelector(selector);
  const canvas = $('#isogeny-canvas');
  const ctx = canvas.getContext('2d');
  const scope = $('#graph-scope');
  const criterion = $('#graph-criterion');
  const search = $('#graph-search');
  const selection = $('#graph-selection');
  const matches = $('#graph-matches');
  let data, nodes, classes, edges, records, openRecord;
  let visible = [], visibleIds = new Set(), visibleEdges = [];
  let scale = 1, offsetX = 0, offsetY = 0, width = 0, height = 0;
  let selectedId = null, dragging = null;

  function layout() {
    let x=30, y=30, rowHeight=0;
    const ordered=[...data.classes].sort((a,b)=>b.size-a.size||a.id.localeCompare(b.id));
    for(const group of ordered) {
      const radius=Math.max(22, 16*Math.sqrt(group.size));
      const diameter=2*radius+46;
      if(x+diameter>1900 && x>30) {x=30;y+=rowHeight;rowHeight=0;}
      group.x=x+radius+23;group.y=y+radius+23;group.radius=radius;
      const ids=[...group.members].sort((a,b)=>Number(nodes.get(b).atlas)-Number(nodes.get(a).atlas)||a.localeCompare(b));
      ids.forEach((id,i)=>{
        const node=nodes.get(id);
        if(ids.length===1) {node.x=group.x;node.y=group.y;return;}
        const angle=i*2.399963229728653;
        const radial=i===0?0:radius*Math.sqrt(i/(ids.length-1));
        node.x=group.x+radial*Math.cos(angle);
        node.y=group.y+radial*Math.sin(angle);
      });
      x+=diameter;rowHeight=Math.max(rowHeight,diameter);
    }
  }

  function refresh() {
    const c=criterion.value;
    const linked=new Set(edges.flatMap(e=>[e.source,e.target]));
    visible=data.nodes.filter(n => (c==='all'||n.checks[c]) &&
      (scope.value==='all'||(scope.value==='explicit'?linked.has(n.id):n.class_size>1||c!=='all')));
    visibleIds=new Set(visible.map(n=>n.id));
    visibleEdges=edges.filter(e=>visibleIds.has(e.source)&&visibleIds.has(e.target));
    if(selectedId&&!visibleIds.has(selectedId)) selectedId=null;
    const classCount=new Set(visible.map(n=>n.class_id)).size;
    $('#graph-stats').textContent=`${visible.length.toLocaleString()} visible curves · ${classCount.toLocaleString()} isogeny classes · ${visibleEdges.length.toLocaleString()} certified directed maps · ${data.counts.graph_nodes.toLocaleString()} total indexed models`;
    fit();showSelection();renderMatches();
  }

  function fit() {
    if(!visible.length) {draw();return;}
    const xs=visible.map(n=>n.x),ys=visible.map(n=>n.y);
    const minX=Math.min(...xs)-45,maxX=Math.max(...xs)+45;
    const minY=Math.min(...ys)-45,maxY=Math.max(...ys)+45;
    scale=Math.min(2.2,Math.max(.08,Math.min(width/(maxX-minX),height/(maxY-minY))*.92));
    offsetX=width/2-(minX+maxX)*scale/2;
    offsetY=height/2-(minY+maxY)*scale/2;
    draw();
  }

  function focus(id) {
    if(!data)return;
    let node=nodes.get(id);
    if(!node)return;
    if(!visibleIds.has(id)) {scope.value='all';criterion.value='all';refresh();node=nodes.get(id);}
    selectedId=id;
    scale=Math.min(2.4,Math.max(.9,scale));
    offsetX=width*.48-node.x*scale;
    offsetY=height*.49-node.y*scale;
    draw();showSelection();
  }

  function toScreen(node) {return {x:node.x*scale+offsetX,y:node.y*scale+offsetY};}
  function draw() {
    ctx.clearRect(0,0,width,height);
    ctx.fillStyle='#0a1725';ctx.fillRect(0,0,width,height);
    const byClass=new Map();
    for(const node of visible) {
      if(!byClass.has(node.class_id))byClass.set(node.class_id,[]);
      byClass.get(node.class_id).push(node);
    }
    for(const [id,members] of byClass) {
      if(members.length<2)continue;
      const group=classes.get(id),center={x:group.x*scale+offsetX,y:group.y*scale+offsetY};
      if(center.x+group.radius*scale<0||center.x-group.radius*scale>width||center.y+group.radius*scale<0||center.y-group.radius*scale>height)continue;
      ctx.strokeStyle='#4a708142';ctx.lineWidth=1;
      ctx.beginPath();ctx.arc(center.x,center.y,Math.max(8,group.radius*scale+8),0,Math.PI*2);ctx.stroke();
      if(scale>.32) {
        ctx.beginPath();for(const node of members) {const p=toScreen(node);ctx.moveTo(center.x,center.y);ctx.lineTo(p.x,p.y);}ctx.stroke();
      }
    }
    ctx.lineWidth=Math.max(1,Math.min(2.5,1.4*scale));
    for(const edge of visibleEdges) {
      const a=toScreen(nodes.get(edge.source)),b=toScreen(nodes.get(edge.target));
      if((a.x<0&&b.x<0)||(a.x>width&&b.x>width)||(a.y<0&&b.y<0)||(a.y>height&&b.y>height))continue;
      ctx.strokeStyle=edge.evidence==='verified_rational_map'?'#edbd83ba':'#75d6e9ab';
      ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke();
    }
    for(const node of visible) {
      const p=toScreen(node);if(p.x<-10||p.y<-10||p.x>width+10||p.y>height+10)continue;
      const radius=node.id===selectedId?6:Math.max(2.3,Math.min(4.2,2.4*scale));
      ctx.fillStyle=node.id===selectedId?'#fff4b8':node.checks.isogeny_transfer?'#edbd83':node.atlas?'#a8f1cf':'#84afcc';
      ctx.beginPath();ctx.arc(p.x,p.y,radius,0,Math.PI*2);ctx.fill();
      if(node.id===selectedId||scale>1.2&&(node.label==='P-192'||node.label==='P-224')) {
        ctx.fillStyle='#f1f7f4';ctx.font='12px ui-sans-serif,system-ui';ctx.fillText(node.label,p.x+8,p.y-9);
      }
    }
    if(!visible.length) {ctx.fillStyle='#a8bdcd';ctx.font='14px ui-sans-serif,system-ui';ctx.fillText('No indexed curves match this criterion and display scope.',24,36);}
  }

  function resize() {
    const rect=canvas.getBoundingClientRect();
    width=Math.max(1,rect.width);height=Math.max(1,rect.height);
    const dpr=Math.min(window.devicePixelRatio||1,2);
    canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);
    ctx.setTransform(dpr,0,0,dpr,0,0);
    fit();
  }

  function nodeAt(clientX,clientY) {
    const rect=canvas.getBoundingClientRect(),x=clientX-rect.left,y=clientY-rect.top;
    let best=null,bestD=12*12;
    for(const node of visible) {
      const p=toScreen(node),d=(p.x-x)**2+(p.y-y)**2;
      if(d<bestD){best=node;bestD=d;}
    }
    return best;
  }

  function appendLine(parent,tag,text,cls) {const el=document.createElement(tag);el.textContent=text;if(cls)el.className=cls;parent.append(el);return el;}
  function showSelection() {
    selection.replaceChildren();
    const node=nodes.get(selectedId);
    if(!node) {appendLine(selection,'p','Select a node or search for a curve to inspect its class and recorded maps.');return;}
    appendLine(selection,'p',node.atlas?'ATLAS MODEL':'VERIFIED WALK MODEL','selection-kicker');
    appendLine(selection,'h3',node.label);
    appendLine(selection,'code',node.id,'selection-id');
    appendLine(selection,'p',`F(${node.p}${node.degree===1?'':`^${node.degree}`}) · trace ${node.trace} · ${node.class_size} indexed model${node.class_size===1?'':'s'} in class`);
    const checks=Object.entries(node.checks);
    if(checks.length) {
      appendLine(selection,'h4','Recorded checks');
      for(const [id,check] of checks) appendLine(selection,'p',`${id} · ${check.status}: ${check.detail}`,'selection-check');
    }
    const incident=edges.filter(e=>e.source===node.id||e.target===node.id);
    appendLine(selection,'h4',`Certified maps (${incident.length})`);
    if(!incident.length) appendLine(selection,'p','No explicit map indexed for this model. Class peers below are related by the trace theorem.');
    for(const edge of incident) {
      const outbound=edge.source===node.id;
      const other=nodes.get(outbound?edge.target:edge.source);
      const row=document.createElement('div');row.className='neighbor-row';
      const button=document.createElement('button');button.type='button';button.className='neighbor-button';
      button.textContent=`degree ${edge.degree} ${outbound?'→':'←'} ${other.label}`;
      button.title=`${edge.evidence} · ${edge.id}`;button.addEventListener('click',()=>focus(other.id));row.append(button);
      const proof=document.createElement('a');proof.href=`https://github.com/aburan28/crypto/blob/codex/weak-curves-leaderboard-20261012/${edge.map_ref||edge.source_ref}`;
      proof.target='_blank';proof.rel='noopener noreferrer';proof.textContent='map ↗';proof.title=edge.certificate_sha256||edge.id;row.append(proof);selection.append(row);
    }
    const peerIds=classes.get(node.class_id).members.filter(id=>id!==node.id&&!incident.some(e=>e.source===id||e.target===id));
    if(peerIds.length) {
      appendLine(selection,'h4',`Other class members (${peerIds.length})`);
      const list=document.createElement('div');list.className='peer-list';selection.append(list);
      for(const id of peerIds) {
        const peer=nodes.get(id),button=document.createElement('button');button.type='button';button.className='peer-button';
        button.textContent=peer.label;button.title='Same field size and trace; explicit map not indexed';
        button.addEventListener('click',()=>focus(id));list.append(button);
      }
    }
    if(node.atlas) {
      const record=records.get(node.id);
      if(record){const button=document.createElement('button');button.type='button';button.className='open-atlas';button.textContent='Open full atlas record ↗';button.addEventListener('click',()=>openRecord(record));selection.append(button);}
    }
    if(incident.length) {
      const refs=[...new Set(incident.map(e=>e.source_ref))];
      for(const ref of refs){const link=document.createElement('a');link.href=`https://github.com/aburan28/crypto/blob/codex/weak-curves-leaderboard-20261012/${ref}`;link.target='_blank';link.rel='noopener noreferrer';link.textContent=`Map evidence: ${ref.split('/').pop()} ↗`;selection.append(link);}
    }
  }

  function renderMatches() {
    if(!data)return;
    matches.replaceChildren();const q=search.value.trim().toLowerCase();
    const found=q?data.nodes.filter(n=>`${n.label} ${n.slug} ${n.id} ${n.class_id}`.toLowerCase().includes(q)).slice(0,12):data.nodes.filter(n=>['P-192','P-224','Oakley Group 3','Oakley Group 4'].includes(n.label));
    if(q&&!found.length){appendLine(matches,'p','No indexed model matches. Try a name, field, or ICV1 fragment.');return;}
    for(const node of found){const button=document.createElement('button');button.type='button';button.textContent=node.label;const small=document.createElement('small');small.textContent=`${node.bits} bits · ${node.class_size} in class`;button.append(small);button.addEventListener('click',()=>focus(node.id));matches.append(button);}
  }

  function zoom(factor,x=width/2,y=height/2) {
    const next=Math.max(.08,Math.min(8,scale*factor));
    offsetX=x-(x-offsetX)*next/scale;offsetY=y-(y-offsetY)*next/scale;scale=next;draw();
  }

  window.addEventListener('weakcurves-ready',event=>{
    data=event.detail.graph;nodes=new Map(data.nodes.map(n=>[n.id,n]));
    classes=new Map(data.classes.map(c=>[c.id,c]));edges=data.edges;
    records=new Map(event.detail.records.map(r=>[r.id,r]));openRecord=event.detail.openRecord;
    for(const item of event.detail.criteria){const option=document.createElement('option');option.value=item.id;option.textContent=`${item.name} (${data.counts.atlas_checks[item.id]||0} atlas)`;criterion.append(option);}
    layout();resize();refresh();
    new ResizeObserver(resize).observe(canvas);
  });
  window.addEventListener('weakcurves-focus-graph',event=>focus(event.detail.id));
  [scope,criterion].forEach(el=>el.addEventListener('change',()=>{if(data)refresh();}));
  search.addEventListener('input',renderMatches);
  $('#graph-overview').addEventListener('click',fit);
  $('#graph-zoom-in').addEventListener('click',()=>zoom(1.4));
  $('#graph-zoom-out').addEventListener('click',()=>zoom(1/1.4));
  canvas.addEventListener('pointerdown',event=>{dragging={x:event.clientX,y:event.clientY,originX:offsetX,originY:offsetY,moved:false};canvas.setPointerCapture(event.pointerId);});
  canvas.addEventListener('pointermove',event=>{if(!dragging)return;const dx=event.clientX-dragging.x,dy=event.clientY-dragging.y;if(Math.abs(dx)+Math.abs(dy)>4)dragging.moved=true;offsetX=dragging.originX+dx;offsetY=dragging.originY+dy;draw();});
  canvas.addEventListener('pointerup',event=>{if(!dragging)return;if(!dragging.moved){const node=nodeAt(event.clientX,event.clientY);if(node){selectedId=node.id;draw();showSelection();}}dragging=null;});
  canvas.addEventListener('wheel',event=>{event.preventDefault();const rect=canvas.getBoundingClientRect();zoom(event.deltaY<0?1.15:1/1.15,event.clientX-rect.left,event.clientY-rect.top);},{passive:false});
})();
