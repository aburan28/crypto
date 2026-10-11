use serde_json::Value;
use std::{fs, path::Path};
fn text(v: &Value, k: &str) -> String {
    v[k].as_str().unwrap_or("").to_string()
}
fn tex(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '\\' => "\\textbackslash{}".into(),
            '{' => "\\{".into(),
            '}' => "\\}".into(),
            '_' => "\\_".into(),
            '#' => "\\#".into(),
            '%' => "\\%".into(),
            '&' => "\\&".into(),
            '$' => "\\$".into(),
            '^' => "\\textasciicircum{}".into(),
            '~' => "\\textasciitilde{}".into(),
            _ => c.to_string(),
        })
        .collect()
}
fn title(s: &str) -> String {
    s.replace('_', " ")
}
fn properties(f: &Value) -> String {
    f["properties"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p.as_str().unwrap())
        .collect::<Vec<_>>()
        .join(" ")
}
fn html(catalogue: &Value, index: &Value) -> String {
    let data = serde_json::to_string(&serde_json::json!({"catalogue":catalogue,"index":index}))
        .unwrap()
        .replace('<', "\\u003c");
    HTML.replace("__DATA__", &data)
}
pub fn check(repo: &Path, catalogue: &Value, index: &Value) {
    assert_eq!(
        fs::read_to_string(repo.join("docs/curves/weak-families/index.html")).unwrap(),
        html(catalogue, index),
        "browser index must match the canonical family and model data"
    );
}
pub fn render(repo: &Path, catalogue: &Value, index: &Value) {
    let out = repo.join("research/weak_curve_catalog_20261010");
    let docs = repo.join("docs/curves/weak-families");
    let families = catalogue["families"].as_array().unwrap();
    let n = index["summary"]["models"].as_u64().unwrap();
    let mut md = format!("# Weak-curve family catalogue\n\nThe first catalogue contains **{} family or condition entries** and **{n} exact model records**. It separates subgroup reductions, covering geometry, isogeny-class support, structural leads and input contracts. Each entry supplies its domain, recognition condition, finite construction/enumeration description, subgroup requirements, source and evidence status.\n\n[Searchable catalogue](../../docs/curves/weak-families/index.html), [canonical family data](../../docs/curves/weak-families/catalog.json), and [model index](../../docs/curves/weak-families/models.json) are generated with the native tool.\n\n![Classification and evidence levels](classification.svg)\n\n", families.len());
    md.push_str("## What is classified\n\nA family-level theorem and an individual model record answer different questions. The catalogue asks whether the stated model condition holds, whether its isogeny class has a family representative, whether that representative can be constructed from the supplied source, and whether a specified subgroup has an established lower-cost route. A result at one level is not promoted to the next.\n\nThe inventories retain their input laws: 277 p7 models come from postselected route, conversion and control replays; 25 larger-degree models come from constructed controls; 512 large models are the original independent source draws. Inventory counts are not family-density estimates. Exact p7 labels use a complete parameter image and a complete ordinary trace oracle with t = 2 modulo 4. The larger source labels remain unresolved after the bounded searches.\n\n");
    md.push_str("## Family records\n\n");
    md.push_str("The published reductions and covering families retain their original attribution. Momose–Chao supply the Type I/II covering classification and its hyperelliptic subcases; Joux–Vitse supply cover/decomposition index calculus. The ISO-1 contribution attached here is the exact ordinary isogeny-class support criterion for the two hyperelliptic branches, the complete finite class oracles, and independently replayed source constructions. The catalogue adds a shared identity and evidence schema across these mechanisms.\n\n");
    for f in families {
        md.push_str(&format!("### {}\n\n`{}` — {}; {}.\n\n**Domain:** {}\n\n**Recognition:** {}\n\n**Find members:** {}\n\n**Completeness:** {}\n\n**Subgroup gate:** {}\n\n**Cost gate:** {}\n\n**Local status:** {}. Sources: {}.\n\n",text(f,"name"),text(f,"id"),title(&text(f,"category")),title(&text(f,"evidence")),text(f,"domain"),text(f,"recognition"),text(f,"find_members"),text(f,"exhaustive_scope"),text(f,"subgroup_gate"),text(f,"cost_gate"),title(&text(f,"local_status")),f["sources"].as_array().unwrap().iter().map(|v| v.as_str().unwrap()).collect::<Vec<_>>().join(", ")));
        md.push_str(&format!("**Properties:** {}\n\n", properties(f)));
    }
    md.push_str("## Inventory counts\n\nThese are counts of exact model records under the stated inventory laws. Model membership and class support counts overlap and must not be added as disjoint populations.\n\n| Label | Records |\n| --- | ---: |\n");
    for (key, value) in index["summary"]["inventory_counts"].as_object().unwrap() {
        md.push_str(&format!("| `{key}` | {value} |\n"));
    }
    md.push_str("\n## Sources\n\n");
    md.push_str("The browser also exposes the complete ordinary p7 trace table separately from the model inventory: 294 eligible classes, 126 cubic-support classes, 210 quadratic-support classes, 88 in both, 248 in the union and 46 zero for these two branches. A zero for this union does not classify other covering families.\n\n");
    for s in catalogue["sources"].as_array().unwrap() {
        let link = s["url"]
            .as_str()
            .map(str::to_string)
            .unwrap_or_else(|| format!("../../{}", s["path"].as_str().unwrap()));
        md.push_str(&format!(
            "- **{}:** [{}]({link}).\n",
            text(s, "id"),
            text(s, "title")
        ));
    }
    md.push_str("\n## Finding all members, and extending the catalogue\n\nFor a fixed field and a parametrized family, exhaustive coverage means enumerating its complete admissible parameter set and quotienting only by documented equivalences. It does not follow from a successful bounded walk. For class support, a complete oracle or a verified representative provides the appropriate label; restricted closures and caps remain unknown. The source and subgroup requirements make the claimed object precise.\n\nThe next catalogue rungs are independent Type I/II recognizers and witnesses, binary-descent fixtures, and trace/order-based baseline fixtures with specified subgroups. Efficient on-demand CM support and verified large-field source routes remain open. Canonical curve IDs and source hashes allow new families to be attached without rewriting earlier measurements.\n\nThe existing IC ratio graphs and two-branch census figures were checked for affected claims. Their numerical values are unchanged: this round organizes established evidence rather than providing a new paired cost. The new classification diagram and inventory table expose the additional catalogue scope.\n\n[Protocol](PROTOCOL.md), [native builder](catalog.rs), and [derived summary](summary.json) retain the input binding and status contract.\n");
    md.push_str("\n## Reproduction and validation\n\nFrom the repository root, use the focused native workspace:\n\n```sh\ncargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bins --lib\ncargo run --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bin weak_curve_catalog -- check .\n```\n\nUse `build .` to regenerate the model and class index, followed by `render .` to regenerate the HTML, Markdown, LaTeX and SVG presentation. The input bindings inside `models.json` freeze five parent evidence files. Checks reject input changes, inconsistent model/trace/order records, duplicate model identities, unsafe evidence promotion and stale browser data. [Validation](VALIDATION.md) records the executed native, browser and PDF checks.\n");
    fs::write(out.join("REPORT.md"), &md).unwrap();
    let mut latex = String::from("\\documentclass[10pt]{article}\n\\usepackage[margin=0.78in]{geometry}\n\\usepackage{amsmath,amssymb,longtable,booktabs,graphicx,xcolor}\n\\usepackage{xurl,needspace}\n\\usepackage[colorlinks=true,breaklinks=true,urlcolor=blue]{hyperref}\n\\emergencystretch=2em\n\\setlength{\\parindent}{0pt}\n\\setlength{\\parskip}{5pt}\n\\title{Weak elliptic-curve families: classification,\\\\ recognition and evidence}\n\\author{Research catalogue and model index}\n\\date{10 October 2026}\n\\begin{document}\n\\maketitle\n");
    latex.push_str(&format!("The catalogue contains {} family or condition entries and {n} exact synthetic model records. It organizes published subgroup reductions and covering families alongside proved class support, verified geometry and structural leads. Every record states its field, subgroup, recognition and construction scope.\n\n",families.len()));
    latex.push_str("\\section{Classification contract}\nA model condition, class-support theorem, source construction and subgroup cost comparison are distinct evidence levels. A covering witness establishes geometry; a subgroup condition and a verified cost are needed for a selected-group comparison. Necessary filters and structural signals retain their own labels.\n\nThe model index contains 277 postselected $p=7$ route/conversion/control models, 25 larger-degree constructed-control models, and 512 original independent large-field sources. The complete $p=7$ torus images and ordinary $t\\equiv2\\pmod4$ trace oracle provide exact local labels. Large-field class labels remain unresolved after the recorded bounded searches. Inventory counts are not population density estimates.\n\n\\begin{figure}[ht]\\centering\\includegraphics[width=\\linewidth]{classification.png}\\caption{Classification proceeds from a scoped property to its appropriate witness and cost contract. Conditional links retain their hypotheses.}\\end{figure}\n\\section{Catalogue entries}\n");
    latex.push_str("The published families retain their original attribution. Momose--Chao classify Type I/II coverings and their hyperelliptic subcases; Joux--Vitse develop cover/decomposition index calculus. The attached ISO-1 contribution is the exact ordinary isogeny-class support criterion for the two hyperelliptic branches, complete finite class oracles, and independently replayed source constructions. This catalogue unifies the identity and evidence records across mechanisms.\n\n");
    for f in families {
        latex.push_str(&format!("\\Needspace{{8\\baselineskip}}\\subsection{{{}}}\n\\texttt{{{}}}\\quad {} / {}.\\par\n\\textbf{{Domain.}} {}\\par\n\\textbf{{Recognition.}} {}\\par\n\\textbf{{Find members.}} {}\\par\n\\textbf{{Completeness.}} {}\\par\n\\textbf{{Subgroup and cost.}} {} {}\\par\n\\textbf{{Local evidence.}} {}. Sources: {}.\n",tex(&text(f,"name")),tex(&text(f,"id")),tex(&title(&text(f,"category"))),tex(&title(&text(f,"evidence"))),tex(&text(f,"domain")),tex(&text(f,"recognition")),tex(&text(f,"find_members")),tex(&text(f,"exhaustive_scope")),tex(&text(f,"subgroup_gate")),tex(&text(f,"cost_gate")),tex(&title(&text(f,"local_status"))),tex(&f["sources"].as_array().unwrap().iter().map(|v|v.as_str().unwrap()).collect::<Vec<_>>().join(", "))));
        latex.push_str(&format!(
            "\\textbf{{Properties.}} {}\\par\n",
            tex(&properties(f))
        ));
    }
    latex.push_str("\\section{Inventory and next rungs}\n\\begin{longtable}{p{0.77\\linewidth}r}\\toprule Label & Records\\\\\\midrule\\endhead\n");
    for (key, value) in index["summary"]["inventory_counts"].as_object().unwrap() {
        latex.push_str(&format!("{} & {}\\\\\n", tex(&title(key)), value));
    }
    latex.push_str("\\bottomrule\\end{longtable}\nThe separate complete p7 class oracle contains 294 eligible traces, 126 cubic-support classes, 210 quadratic-support classes, 88 in both, 248 in the union, and 46 zeros for this two-branch union.\n\nThese counts have overlapping predicates and separate input laws. The next implementations are independently validated Type I/II recognizers, binary-descent and trace/order fixtures with specified subgroups, and on-demand large-class support. The existing measured ratio graphs retain their numerical values. This round adds classification, source links and exact instance organization.\n\\section{Sources}\n");
    for s in catalogue["sources"].as_array().unwrap() {
        let link = s["url"].as_str().map(str::to_string).unwrap_or_else(|| {
            format!(
                "https://github.com/aburan28/crypto/blob/codex/iso1-report-final-20261008/{}",
                s["path"].as_str().unwrap()
            )
        });
        latex.push_str(&format!(
            "\\Needspace{{3\\baselineskip}}\\textbf{{{}}}: \\href{{{}}}{{{}}}.\\par\n",
            tex(&text(s, "id")),
            link,
            tex(&text(s, "title"))
        ));
    }
    latex.push_str("\\end{document}\n");
    fs::write(out.join("REPORT.tex"), latex).unwrap();
    fs::write(out.join("classification.svg"), DIAGRAM).unwrap();
    fs::write(docs.join("index.html"), html(catalogue, index)).unwrap();
    fs::write(docs.join("README.md"),format!("# Weak-curve families\n\n[Browse the catalogue](index.html): {} family or condition entries and {n} exact synthetic model records. The [family JSON](catalog.json) records hypotheses, recognition, construction, sources and evidence status. The complete 294-class p7 oracle is searchable alongside the models. The [model index](models.json) preserves exact ICV1 identities, model parameters, source bindings, model/class labels and unknown subgroup costs.\n\n[Report and reproduction](../../../research/weak_curve_catalog_20261010/REPORT.md) and [protocol](../../../research/weak_curve_catalog_20261010/PROTOCOL.md) describe the evidence contract. Inventory counts are not empirical family densities. Exact p7 labels, larger-degree geometry controls and unresolved independent large-field sources are separate scopes.\n\nAdd a family with its field/subgroup domain, necessary or sufficient predicates, a construction or finite enumeration description, completeness boundary, source and validation status. Add an instance under its existing ICV1/EC1 identity with immutable evidence; assign each property independently.\n",families.len())).unwrap();
}

const DIAGRAM: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="860" viewBox="0 0 1400 860">
<rect width="1400" height="860" fill="white"/><style>text{font-family:Arial,sans-serif;fill:#17304c}.head{font-size:29px;font-weight:bold}.label{font-size:23px;font-weight:bold}.body{font-size:20px}.small{font-size:18px}rect.box{stroke:#54759a;stroke-width:2;rx:12}</style>
<text class="head" x="55" y="54">Weak-curve catalogue: properties, witnesses and discovery scope</text>
<text class="body" x="55" y="91">21 family/condition entries · 814 exact models · model membership and class support are separate labels</text>
<rect class="box" x="55" y="130" width="410" height="205" fill="#ecf3ff"/>
<text class="label" x="76" y="167">Group and transfer conditions</text><text class="body" x="76" y="205">Anomalous, selected order, pairings</text><text class="body" x="76" y="240">Supersingularity, actual subfield points</text><text class="small" x="76" y="281">Require the specified subgroup and field.</text><text class="small" x="76" y="308">A structural reduction has a cost contract.</text>
<rect class="box" x="495" y="130" width="410" height="205" fill="#e6f6ef"/>
<text class="label" x="516" y="167">Covering families</text><text class="body" x="516" y="205">Cubic / quadratic genus-3 branches</text><text class="body" x="516" y="240">Type I / II, binary and odd-degree lifts</text><text class="small" x="516" y="281">Parameter domain → verified model / cover</text><text class="small" x="516" y="308">Genus and subgroup transfer stay explicit.</text>
<rect class="box" x="935" y="130" width="410" height="205" fill="#fff3df"/>
<text class="label" x="956" y="167">Leads and input contracts</text><text class="body" x="956" y="205">CM, conductor, endomorphism, torsion</text><text class="body" x="956" y="240">Twist validation and model validity</text><text class="small" x="956" y="281">A lead needs a family-specific witness.</text><text class="small" x="956" y="308">Protocol conditions keep their own scope.</text>
<path d="M700 335 V380" stroke="#54759a" stroke-width="3"/><path d="M685 367 L700 384 L715 367" fill="none" stroke="#54759a" stroke-width="3"/>
<rect class="box" x="260" y="390" width="880" height="137" fill="#eff2f8"/><text class="label" x="284" y="428">Isogeny-class support and construction from an independent source</text><text class="body" x="284" y="466">Exact family oracle or explicit representative → class label; verified route → source construction</text><text class="small" x="284" y="501">Restricted closure or finite cap → unresolved. Verify subgroup map and destination cost separately.</text>
<path d="M700 527 V574" stroke="#54759a" stroke-width="3"/><path d="M685 562 L700 578 L715 562" fill="none" stroke="#54759a" stroke-width="3"/>
<rect class="box" x="55" y="588" width="1290" height="155" fill="#fafbfc"/><text class="label" x="80" y="627">Concrete index: exact model IDs, trace/order, sources, and independent property labels</text><text class="body" x="80" y="670">277 p7 route/conversion/control models + 25 larger-degree controls + 512 independent large sources</text><text class="body" x="80" y="709">Complete p7 model/class oracles · geometric controls · unresolved large-field support · subgroup costs unmeasured</text>
<text class="small" x="55" y="800">Inventory counts have distinct input laws. Families can overlap; their counts are not summed as disjoint populations.</text><text class="small" x="55" y="830">Sources: catalog.json, models.json, the frozen two-branch census, coordinate replays and independent population receipts.</text>
</svg>"##;

const HTML: &str = r##"<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Weak curve catalogue</title>
<style>:root{color-scheme:light dark;--bg:#f6f8fb;--ink:#17304c;--card:#fff;--line:#dce3ed;--muted:#516379;--accent:#245aaa}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 system-ui,sans-serif}main{max-width:1230px;margin:auto;padding:30px 24px}h1{font-size:35px;margin:0 0 8px}h2{font-size:22px}.intro{max-width:880px;color:var(--muted)}.bar{display:flex;gap:10px;flex-wrap:wrap;margin:22px 0}input,select,button{font:inherit;border:1px solid var(--line);border-radius:7px;padding:9px 12px;background:var(--card);color:var(--ink)}input{flex:1;min-width:200px}button.active{background:var(--accent);color:white}.cards{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:19px}.card h2{margin:0 0 10px}.badge{display:inline-block;font-size:13px;color:var(--muted);background:var(--bg);padding:3px 8px;border-radius:5px;margin:0 5px 7px 0}.card p{margin:8px 0}.small{font-size:13px;color:var(--muted)}details{margin-top:12px}summary{cursor:pointer;font-weight:600}a{color:var(--accent)}code{font:12px/1.5 ui-monospace,monospace;word-break:break-all}.scroll{overflow:auto;background:var(--card);border:1px solid var(--line);border-radius:8px}table{border-collapse:collapse;width:100%;min-width:740px}td,th{padding:11px 14px;text-align:left;vertical-align:top;border-bottom:1px solid var(--line)}td:first-child{max-width:370px;overflow-wrap:anywhere}.hide{display:none}.note{padding:13px 16px;border-left:3px solid #7092ba;background:var(--card)}.detail{background:var(--card);padding:20px;border:1px solid var(--line);border-radius:9px;margin:18px 0;overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}@media(max-width:700px){main{padding:20px 14px}h1{font-size:28px}.cards{grid-template-columns:1fr}.bar>*{max-width:100%}input{min-width:100%;width:100%}}@media(prefers-color-scheme:dark){:root{--bg:#111b29;--ink:#e4edf8;--card:#182638;--line:#34455c;--muted:#b1bfd3;--accent:#8cb9ff}}</style>
<main><h1>Weak curve catalogue</h1><p class="intro">Families, recognition conditions and construction methods, with exact model records and evidence. Model membership, isogeny-class support and selected-subgroup cost have separate labels.</p><p id="stats" class="small"></p><div class="bar"><button id="familiesTab" class="active">Families</button><button id="curvesTab">Curves</button><button id="classesTab">Classes</button><input id="search" type="search" placeholder="Search a family, property, field or curve ID" aria-label="Search catalogue"><select id="filter" aria-label="Catalogue filter"></select></div><p id="count" class="small"></p><div id="families" class="cards"></div><section id="curves" class="hide"><p class="note">The replay/control inventory and independent large sources have different input laws. These counts are an inventory, rather than a family-density estimate. Unclassified and unresolved labels retain their meaning.</p><div id="curveDetail" class="hide detail"></div><div class="scroll"><table><thead><tr><th>Exact model</th><th>Field bits / degree</th><th>Inventory</th><th>Class evidence</th><th>Direct genus-3 model</th></tr></thead><tbody id="rows"></tbody></table></div><div class="bar"><button id="previous">Previous</button><span id="page"></span><button id="next">Next</button></div></section><section id="classes" class="hide"><p class="note">Complete ordinary p7 oracle: 294 traces with t = 2 modulo 4; 248 have support in the two recorded genus-3 branches and 46 have zero support in that union. The table classifies these two branches only.</p><div class="scroll"><table><thead><tr><th>Field</th><th>Trace</th><th>Order</th><th>Cubic support</th><th>Quadratic support</th><th>Combined support</th></tr></thead><tbody id="classRows"></tbody></table></div></section><details><summary>Sources, scope and growth</summary><p>Coverage is exact only in each declared family and domain. Family overlaps are retained. A geometric witness, structural signal and verified subgroup cost answer different questions.</p><ul id="sources"></ul><ul id="growth"></ul><p><a href="catalog.json">Family JSON</a> · <a href="models.json">Exact model index</a> · <a href="../../../research/weak_curve_catalog_20261010/REPORT.md">Report and reproduction</a></p></details><noscript><p>The catalogue JSON, exact model index and report are available in the source links above. Filtering requires JavaScript.</p></noscript></main>
<script id="data" type="application/json">__DATA__</script><script>
const data=JSON.parse(document.getElementById('data').textContent),catalogue=data.catalogue,index=data.index,$=s=>document.getElementById(s);let tab='families',page=0,filtered=[];const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));const pretty=s=>String(s??'').replaceAll('_',' ');const sources=new Map(catalogue.sources.map(s=>[s.id,s]));const link=s=>s.url||'../../../'+s.path;const readable=s=>({'yes':'Recognized','no':'Absent in complete local image','unresolved':'Unresolved','geometric_support_exact':'Exact geometric support','exact_family_zero':'Exact zero for this family','outside_domain':'Outside recorded genus-3 domain','excluded_by_proved_branch_condition':'Excluded by cubic condition','unclassified_in_this_index':'Unclassified','unclassified':'Unclassified','geometric_control_witness':'Constructed geometry control','unmeasured':'Unmeasured'}[s]||pretty(s));
$('stats').textContent=`${catalogue.families.length} family/condition entries · ${index.models.length} exact models · ${index.class_oracle.length} exact class records · 512 independent large-field sources`;
function filters(){const values=tab==='families'?[...new Set(catalogue.families.map(f=>f.category))]:tab==='classes'?['supported','zero','cubic','quadratic','both']:['independent_model_law_large_source','postselected_replayed_route_or_control','unresolved','geometric_support_exact','exact_family_zero'];$('filter').innerHTML='<option value="">All</option>'+values.map(v=>`<option value="${esc(v)}">${esc(pretty(v))}</option>`).join('');}
function draw(){const q=$('search').value.toLowerCase(),filter=$('filter').value;if(tab==='classes'){drawClasses(q,filter);return;}if(tab==='families'){const all=catalogue.families.filter(f=>(!filter||f.category===filter)&&JSON.stringify(f).toLowerCase().includes(q));$('count').textContent=`${all.length} family/condition entries`;$('families').innerHTML=all.map(f=>`<article class="card"><h2>${esc(f.name)}</h2><span class="badge">${esc(pretty(f.category))}</span><span class="badge">${esc(pretty(f.evidence))}</span><p>${esc(f.domain)}</p><p><strong>Recognition:</strong> ${esc(f.recognition)}</p><details><summary>Find members and verify scope</summary><p><strong>Find members:</strong> ${esc(f.find_members)}</p><p><strong>Completeness:</strong> ${esc(f.exhaustive_scope)}</p><p><strong>Subgroup:</strong> ${esc(f.subgroup_gate)}</p><p><strong>Cost:</strong> ${esc(f.cost_gate)}</p><ul>${f.properties.map(p=>`<li>${esc(p)}</li>`).join('')}</ul><p class="small">${esc(pretty(f.local_status))} · <code>${esc(f.id)}</code></p><p>${f.sources.map(id=>`<a href="${esc(link(sources.get(id)))}">${esc(id)}</a>`).join(' · ')}</p></details></article>`).join('');}else{filtered=index.models.filter(c=>(!filter||c.inventory===filter||c.labels.combined_class===filter||c.labels.higher_family_class===filter)&&JSON.stringify(c).toLowerCase().includes(q));page=Math.max(0,Math.min(page,Math.max(0,Math.ceil(filtered.length/50)-1)));$('count').textContent=`${filtered.length} exact model records`;$('rows').innerHTML=filtered.slice(page*50,(page+1)*50).map(c=>`<tr><td><button class="model" data-id="${esc(c.icv1)}">${esc(c.slug)}</button></td><td>${c.field.bit_length} / ${c.field.degree}</td><td>${esc(pretty(c.inventory))}</td><td>${esc(readable(c.labels.combined_class==='outside_domain'?c.labels.higher_family_class:c.labels.combined_class))}</td><td>Cubic: ${esc(readable(c.labels.model_g3_cubic))}<br>Quadratic: ${esc(readable(c.labels.model_g3_quadratic))}</td></tr>`).join('');$('page').textContent=`Page ${page+1} / ${Math.max(1,Math.ceil(filtered.length/50))}`;$('previous').disabled=page===0;$('next').disabled=(page+1)*50>=filtered.length;document.querySelectorAll('.model').forEach(b=>b.onclick=()=>detail(index.models.find(c=>c.icv1===b.dataset.id)));}}
function detail(c){$('curveDetail').classList.remove('hide');$('curveDetail').innerHTML=`<h2>${esc(c.slug)}</h2><p><code>${esc(c.icv1)}</code></p><p><strong>Trace:</strong> ${esc(c.trace)}<br><strong>Order:</strong> ${esc(c.order)}<br><strong>Subgroup and comparative cost:</strong> unspecified / ${esc(readable(c.labels.subgroup_advantage))}</p><details open><summary>Exact field and model</summary><pre>${esc(c.model_json)}</pre></details><details><summary>Property labels and evidence</summary><pre>${esc(JSON.stringify({labels:c.labels,evidence:c.evidence},null,2))}</pre></details>`;$('curveDetail').scrollIntoView({block:'nearest'});}
function drawClasses(q,filter){const rows=index.class_oracle.filter(c=>(!filter||(filter==='zero'&&!c.combined_support)||(filter==='supported'&&c.combined_support)||(filter==='cubic'&&c.cubic_support)||(filter==='quadratic'&&c.quadratic_support)||(filter==='both'&&c.cubic_support&&c.quadratic_support))&&JSON.stringify(c).toLowerCase().includes(q));$('count').textContent=`${rows.length} exact class records`;$('classRows').innerHTML=rows.map(c=>`<tr><td>7^6</td><td>${esc(c.trace)}</td><td>${esc(c.order)}</td><td>${c.cubic_support?'Yes':'Zero'}</td><td>${c.quadratic_support?'Yes':'Zero'}</td><td>${c.combined_support?'Yes':'Zero'}</td></tr>`).join('');}
function setTab(t){tab=t;page=0;$('classes').classList.toggle('hide',t!=='classes');$('classesTab').classList.toggle('active',t==='classes');$('families').classList.toggle('hide',t!=='families');$('curves').classList.toggle('hide',t!=='curves');$('familiesTab').classList.toggle('active',t==='families');$('curvesTab').classList.toggle('active',t==='curves');filters();draw();}
$('familiesTab').onclick=()=>setTab('families');$('curvesTab').onclick=()=>setTab('curves');$('classesTab').onclick=()=>setTab('classes');$('search').oninput=()=>{page=0;draw()};$('filter').onchange=()=>{page=0;draw()};$('previous').onclick=()=>{page--;draw()};$('next').onclick=()=>{page++;draw()};$('sources').innerHTML=catalogue.sources.map(s=>`<li><a href="${esc(link(s))}">${esc(s.title)}</a></li>`).join('');$('growth').innerHTML=catalogue.growth_requirements.map(s=>`<li>${esc(s)}</li>`).join('');filters();draw();
</script></html>"##;
