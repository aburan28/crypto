//! Standalone property browser; earlier model/family artefacts stay unchanged.
use serde_json::Value;
use std::{fs, path::Path};
fn html(data: &Value) -> String {
    HTML.replace(
        "__DATA__",
        &serde_json::to_string(data).unwrap().replace('<', "\\u003c"),
    )
}
pub fn render(repo: &Path, data: &Value) {
    fs::write(
        repo.join("docs/curves/weak-families/properties.html"),
        html(data),
    )
    .unwrap();
    report(repo, data);
}
fn count(data: &Value, cohort: &str, key: &str) -> u64 {
    data["summary"]["cohorts"][cohort][key]
        .as_u64()
        .unwrap_or(0)
}
fn exact(data: &Value, cohort: &str) -> u64 {
    data["summary"]["cohorts"][cohort]
        .as_object()
        .unwrap()
        .iter()
        .filter(|(k, _)| k.starts_with("endomorphism_v2:"))
        .map(|(_, v)| v.as_u64().unwrap())
        .sum()
}
fn report(repo: &Path, data: &Value) {
    let out = repo.join("research/weak_curve_properties_20261010");
    let cohorts = [
        ("p7_replay_control", "p7 replay/control", "p7"),
        ("higher_degree_control", "Higher-degree controls", "higher"),
        ("independent_large", "Independent large sources", "large"),
    ];
    let mut table=String::from("| Inventory | Models | Full 2 | Full 4 | Exact 2-level | Bounded 2-level | Proper j | Proper displayed model | Proper isomorphic model |\n| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n");
    let mut tex_table=String::from("\\begin{tabular}{lrrrrrr}\\toprule Inventory & Models & Full 2 & Full 4 & Exact level & Bounded & Proper iso\\\\\\midrule\n");
    let mut tsv = String::from("cohort exact bounded\n");
    for (c, name, tag) in cohorts {
        let n = count(data, c, "models");
        let e = exact(data, c);
        let b = count(data, c, "endomorphism_v2_unknown");
        table.push_str(&format!(
            "| {name} | {n} | {} | {} | {e} | {b} | {} | {} | {} |\n",
            count(data, c, "full_2_torsion"),
            count(data, c, "full_4_torsion"),
            count(data, c, "j_proper_subfield"),
            count(data, c, "literal_proper_subfield"),
            count(data, c, "isomorphic_proper_subfield")
        ));
        tex_table.push_str(&format!(
            "{name} & {n} & {} & {} & {e} & {b} & {}\\\\\n",
            count(data, c, "full_2_torsion"),
            count(data, c, "full_4_torsion"),
            count(data, c, "isomorphic_proper_subfield")
        ));
        tsv.push_str(&format!("{tag} {e} {b}\n"));
    }
    tex_table.push_str("\\bottomrule\\end{tabular}\n");
    fs::write(out.join("level_coverage.tsv"), tsv).unwrap();
    fs::write(out.join("level_coverage.tex"), PLOT).unwrap();
    fs::write(out.join("property_pipeline.svg"), DIAGRAM).unwrap();
    let certificates: Value =
        serde_json::from_slice(&fs::read(out.join("identity_certificates.json")).unwrap()).unwrap();
    let ids = certificates
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c["id"].as_str().unwrap())
        .collect::<Vec<_>>()
        .join(", ");
    let mut md=String::from("# Certified structural properties of the curve catalogue\n\nNative arithmetic and independent PARI agree on **814 exact models**. The catalogue now certifies **360 ordinary 2-volcano levels** and bounds the remaining **454**. The 512 independent 192–252-bit sources include **219 certified 2-surface models** and **65 models with full rational 4-torsion**. Every property record retains its original family/class labels and exact model identity.\n\n[Search the property browser](../../docs/curves/weak-families/properties.html), [canonical property records](../../docs/curves/weak-families/properties.json), [protocol](PROTOCOL.md), and [validation](VALIDATION.md).\n\n![Property verification pipeline](property_pipeline.svg)\n\n## Inventory and field scope\n\nThe p7 inventory contains 277 route/conversion/control models; the 25 higher-degree models are constructed controls. The 512 independent large sources use the original model law, including full rational 2-torsion by construction. All 814 retained records are ordinary. These exact counts describe the frozen inventories and have no uniform-class density or sampling-error interpretation. Field sizes are 7^6 (17 bits), and 192–252 bits in degrees 6, 10 and 14. Selected subgroup orders and generators remain unspecified.\n\n");
    md.push_str(&table);
    md.push_str("\n![Exact versus bounded ordinary 2-levels](level_coverage.svg)\n\n## Certified properties\n\nEach record contains ordinary/supersingular status, exact ambient-order valuations at 2,3,5,7,11,13, the unclassified residual cofactor, Frobenius discriminant, Frobenius conductor 2-depth, splitting of 2 in the CM field, rational nonidentity 2-point count, full rational 4-torsion, three field degrees, a twist power-class witness, and exact or bounded local endomorphism depth. Ten p7 records have a fully 13-smooth ambient order. Larger ambient cofactors are not silently treated as prime, and selected-point containment remains unspecified.\n\nThe property checks supply scoped supersingular, prime-field anomalous, ambient smoothness, model-descent and necessary genus-3 condition labels. A passed necessary condition retains that status. Existing exact p7 family labels and unresolved large-field class labels are unchanged.\n\n## Local conductor and splitting without odd factorization\n\nWrite Delta=t^2-4Q=2^s u with u odd, Delta<0 for an ordinary model. If s is odd, v2(D_K)=3. If s is even and u=1 mod4, v2(D_K)=0; if u=3 mod4, v2(D_K)=2. Thus v2(f_pi)=(s-v2(D_K))/2. Odd squares are 1 mod8, so 2 splits when the resulting odd fundamental-discriminant residue is 1 mod8, is inert at 5 mod8, and is ramified when its 2-valuation is positive. This computes local data while full D_K and the odd conductor remain unknown.\n\nThe independent method searches for integral (pi-a)/2^k. Its trace and norm must be integers: (t-2a)/2^k and (Q-at+a^2)/2^(2k). For even t there are two candidate a residues modulo 2^k; for odd t no positive k works. The maximum k is the local conductor depth, since pi=A+f_pi omega in a maximal-order basis. PARI checks this method on every ordinary source and also checks the full fundamental discriminant for the p7-sized discriminants.\n\n## Proposition: torsion certifies ordinary 2-levels\n\nLet E/F_Q be ordinary in odd characteristic, with End(E)=Z+f_E O_K and s=v2(f_pi). If E has exactly one rational nonidentity 2-point, then v2(f_E)=s: the model is on the ordinary 2-volcano floor. Full rational 2-torsion gives v2(f_E)<=s-1; full rational 4-torsion gives v2(f_E)<=s-2. When a bound collapses to zero, the local level is exact. Full endomorphism conductors remain uncomputed.\n\nProof: choose a basis 1,omega of End(E), and write pi=a+b omega. The order index gives f_pi=f_E |b|. Full E[2] rational is equivalent to pi-1 lying in 2 End(E), because [2] is separable and an endomorphism killing E[2] factors through [2]. If b is even then Norm(pi)=Q odd forces a odd, so pi-1 lies in 2 End(E). Hence absence of full rational 2-torsion forces b odd and gives the floor equality. Full E[4] rational similarly implies pi-1 lies in 4 End(E), giving the second upper bound. The degree-2 volcano model applies because the characteristic exceeds three.\n\nNative arithmetic obtains the rational 2-point count from the ambient order parity and cubic discriminant square class, then validates every PARI-supplied root. Full 4-torsion is certified from the root differences and -1 square class. PARI uses the sextic part of the degree-4 division polynomial and its y-coordinate square conditions independently.\n\n## Three fields of definition\n\nThe minimal Frobenius period of j is its field degree. The joint period of a,b is the displayed model's degree. To test the smallest field admitting a K-isomorphic model, use the finite cyclic twist quotient. For generic j, normalize a0=3j/(1728-j), b0=2j/(1728-j), with twist parameter D=3b/(2a). For j=0 use b and sixth powers; for j=1728 use a and fourth powers.\n\nFor a candidate F_(p^d), put M=p^n-1, m=p^d-1 and h=gcd(twist_order,M/m). A twist parameter lies in the image of the subfield twist classes precisely when its M/h power is one: F_(p^d)^* times (K^*)^twist_order is generated by the h-th power of a cyclic generator. Each record retains all candidate-field power witnesses, including rejected smaller fields. This proves the minimal degree without providing a selected-point subgroup containment claim.\n\nFour p7 records have j degree three. Three admit an isomorphic model over F_(7^3), while one requires F_(7^6). Only one already displays coefficients over F_(7^3).\n\n| Exact model | j degree | Displayed coefficient degree | K-isomorphic-model degree |\n| --- | ---: | ---: | ---: |\n| `icv1-fp3k6-tm110-5f712c89` | 3 | 6 | 3 |\n| `icv1-fp3k6-t202-39e02f24` | 3 | 6 | 6 |\n\n## Certificates and independent replay\n\n");
    md.push_str(&format!("Thirteen portable identity records ({ids}) cover the cubic discriminant, cleared normalized-j and twist coefficients, translated Frobenius norm, and the eight unit-group index degrees used here. Each has 64 saved and 128 fresh evaluations; native expansion also verifies exact integer coefficient equality. Their formal-proof pointer remains null. The order/torsion and cyclic-group propositions have deductive proofs above.\n\n"));
    md.push_str("PARI audit run 2 completes 814 rows with empty stderr under the 300-second cap. Run 1 stopped during stack-setting initialization; its zero-row receipt and executed sources remain archived. Source hashes, wrapper hashes, fixed input hash, runtime metadata and raw output are retained. Point counts are reused from the independently verified parent evidence, rather than claimed as newly recomputed here.\n\nThe previous catalogue study and two-branch evidence remain sealed. The new property browser is a sidecar; old model/family data are unchanged. Existing IC ratio and census graphs retain their values because this round adds structure, not a new paired cost.\n\n## Reproduction\n\n```sh\ncargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bins --lib\ncargo run --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bin weak_curve_properties -- check .\n```\n\nThe `emit` operation recreates fixed GP inputs. `audit` preserves a fresh run directory and refuses an existing one; archived runs are independently replayable through `replay.gp` with `WEAK_PROPERTIES_INPUTS` set to the frozen inputs. `build` rebuilds the property sidecar; `render` regenerates the browser, report and figures.\n\nSources: [Sutherland, ordinary/supersingular curves and endomorphism orders](https://math.mit.edu/classes/18.783/2023/LectureNotes13.pdf), [Sutherland, isogeny volcanoes](https://math.mit.edu/classes/18.783/2023/LectureNotes22.pdf), [PARI arithmetic functions](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Arithmetic_functions.html), [PARI finite-field polynomial roots](https://pari.math.u-bordeaux.fr/dochtml/html/Polynomials_and_power_series.html), and [parent theorem/data](../iso1_weak_classes_20261007/two_branch_20261009/REPORT.md).\n");
    fs::write(out.join("REPORT.md"), md).unwrap();
    fs::write(out.join("REPORT.tex"), TEX.replace("__TABLE__", &tex_table)).unwrap();
}
const PLOT: &str = r#"\documentclass[tikz,border=8pt]{standalone}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}
\begin{document}\begin{tikzpicture}\begin{axis}[ybar stacked,width=14cm,height=7cm,ymin=0,ymax=590,ylabel={Exact model records},symbolic x coords={p7,higher,large},xtick=data,xticklabels={p7 replay/control,Higher-degree controls,Independent large},x tick label style={font=\small},bar width=45pt,legend style={at={(0.5,1.02)},anchor=south,legend columns=2},nodes near coords,point meta=explicit,every node near coord/.append style={font=\footnotesize},grid=major,major grid style={gray!20}]
\addplot[fill=blue!65,draw=blue!65] table[x=cohort,y=exact,meta=exact] {level_coverage.tsv};
\addplot[fill=gray!40,draw=gray!50] table[x=cohort,y=bounded,meta=bounded] {level_coverage.tsv};
\legend{Certified ordinary 2-level,Bounded ordinary 2-level}
\end{axis}\end{tikzpicture}\end{document}
"#;
const DIAGRAM: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="720" viewBox="0 0 1400 720"><rect width="1400" height="720" fill="white"/><style>text{font-family:Arial,sans-serif;fill:#17304c}.title{font-size:29px;font-weight:bold}.head{font-size:24px;font-weight:bold}.body{font-size:21px}.small{font-size:19px}rect.box{stroke:#52779e;stroke-width:2;rx:12}</style>
<text x="50" y="52" class="title">Certified curve properties: local levels, torsion and fields of definition</text><text x="50" y="91" class="body">814 exact models · fixed source hash · original family/class labels and inventory laws retained</text>
<rect class="box" x="50" y="135" width="410" height="185" fill="#eaf2ff"/><text x="73" y="171" class="head">Trace and ordinary CM order</text><text x="73" y="214" class="body">Frobenius discriminant and 2-depth</text><text x="73" y="248" class="body">Splitting of 2, small order factors</text><text x="73" y="290" class="small">Independent trace/norm integrality scan</text>
<rect class="box" x="495" y="135" width="410" height="185" fill="#e8f5ef"/><text x="518" y="171" class="head">Rational torsion and local levels</text><text x="518" y="214" class="body">All rational 2-roots; full 4-torsion</text><text x="518" y="248" class="body">Exact 2-floor or certified level bounds</text><text x="518" y="290" class="small">Root witnesses and division polynomial</text>
<rect class="box" x="940" y="135" width="410" height="185" fill="#fff2dc"/><text x="963" y="171" class="head">Three fields of definition</text><text x="963" y="214" class="body">j / displayed coefficients / isomorphism</text><text x="963" y="248" class="body">Frobenius periods and twist classes</text><text x="963" y="290" class="small">Subfield power witnesses keep each scope</text>
<path d="M700 320 V370" stroke="#52779e" stroke-width="3"/><path d="M688 355 L700 372 L712 355" fill="none" stroke="#52779e" stroke-width="3"/>
<rect class="box" x="190" y="385" width="1020" height="140" fill="#eef1f7"/><text x="213" y="425" class="head">Native arithmetic + independent PARI + identity certificates</text><text x="213" y="467" class="body">360 exact ordinary 2-levels · 454 bounded levels · full endomorphism conductors remain unknown</text><text x="213" y="503" class="small">Ordinary 2-surface is a local statement at 2; it does not assign maximality at other primes.</text>
<rect class="box" x="50" y="568" width="1300" height="107" fill="#fafbfc"/><text x="75" y="610" class="head">Catalogue properties preserve their evidence boundaries</text><text x="75" y="650" class="body">Necessary conditions / family membership / class support / selected-point containment / subgroup cost</text></svg>"##;
const TEX: &str = r#"\documentclass[10pt]{article}
\usepackage[margin=0.75in]{geometry}\usepackage{amsmath,amssymb,booktabs,graphicx,xurl,needspace}
\usepackage[colorlinks=true,urlcolor=blue]{hyperref}\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}\emergencystretch=2em
\title{Certified structural properties of the curve catalogue}\author{Fixed-model classification and independent replay}\date{11 October 2026}
\begin{document}\maketitle
Native arithmetic and PARI agree on all 814 exact models. We certify 360 ordinary 2-volcano levels and retain bounds for 454. Among the 512 independent 192--252-bit sources, 219 are certified at the 2-surface and 65 have full rational 4-torsion. Original model identities and family/class labels are retained.
\section{Inventory and evidence}
The 277 $p=7$ route/conversion/control models are postselected; 25 higher-degree models are constructed controls. The 512 independent large sources retain their original model law, including full rational 2-torsion by construction. All retained models are ordinary. Counts describe these inventories, with no uniform-class density or sampling-error interpretation. Field sizes are $7^6$ (17 bits) and 192--252 bits, in degrees 6, 10 and 14. Selected subgroups and generators remain unspecified.
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{property_pipeline.png}\caption{Independent structural properties feed the existing catalogue without promoting a necessary filter to family membership.}\end{figure}
\begin{center}\small __TABLE__\end{center}
\begin{figure}[ht]\centering\includegraphics[width=0.92\linewidth]{level_coverage.pdf}\caption{Exact and bounded ordinary 2-levels in each frozen model cohort. Units are model records; cohort input laws differ.}\end{figure}
\section{Local conductor data without odd factorization}
For an ordinary trace write $\Delta=t^2-4Q=2^s u<0$ with $u$ odd. Fundamental-discriminant classification gives
\[
v_2(D_K)=\begin{cases}3&s\text{ odd},\\0&s\text{ even and }u\equiv1\pmod4,\\2&s\text{ even and }u\equiv3\pmod4.\end{cases}
\qquad v_2(f_\pi)=\frac{s-v_2(D_K)}2.
\]
Odd squares are one modulo eight, so the odd fundamental-discriminant residue determines split/inert behavior of two; a positive two-valuation is ramified. Full $D_K$ and the odd conductor remain uncomputed.

PARI independently searches integral $(\pi-a)/2^k$, whose trace and norm must both be integers:
\[
\frac{t-2a}{2^k},\qquad \frac{Q-at+a^2}{2^{2k}}.
\]
For even $t$ there are two possible $a$ residues modulo $2^k$; odd $t$ permits no positive $k$. In a maximal-order basis $\pi=A+f_\pi\omega$, the largest admissible denominator is exactly $v_2(f_\pi)$. This method checks every ordinary model, with a full-discriminant cross-check on the small discriminants.
\section{Torsion certifies ordinary 2-levels}
\textbf{Proposition.} Let $E/\mathbb F_Q$ be ordinary in odd characteristic and let $s=v_2(f_\pi)$. Exactly one rational nonidentity 2-point implies $v_2(f_E)=s$, the ordinary 2-volcano floor. Full rational 2-torsion implies $v_2(f_E)\le s-1$; full rational 4-torsion implies $v_2(f_E)\le s-2$. A collapsed bound certifies the local level.

\textbf{Proof.} In a basis $1,\omega$ of $\operatorname{End}(E)$, write $\pi=a+b\omega$. The order index is $f_\pi=f_E|b|$. Full rational $E[2]$ is equivalent to $\pi-1\in2\operatorname{End}(E)$: an endomorphism killing $E[2]$ factors through the separable map $[2]$. If $b$ is even, the odd norm $Q$ forces $a$ odd, giving this divisibility. Absence of full rational 2-torsion therefore forces $b$ odd and the floor equality. Full rational $E[4]$ similarly implies divisibility by four and the second bound.\hfill$\square$

The native checker validates all rational 2-roots and uses their differences and the square class of $-1$ to certify full 4-torsion. PARI independently uses the sextic factor of the degree-four division polynomial and the square classes of its $y$ coordinates. Whole endomorphism conductors remain unknown; a 2-surface statement does not assert maximality at other primes.
\section{Three fields of definition}
The least Frobenius period of $j$ determines its degree; the joint period of $a,b$ determines the displayed model's degree. A $K$-isomorphic model can have a different minimum degree.

For generic $j$, set $a_0=3j/(1728-j)$ and $b_0=2j/(1728-j)$. The twist parameter is $D=3b/(2a)$, with $a=D^2a_0$ and $b=D^3b_0$. At $j=0$ use $b$ modulo sixth powers, and at $j=1728$ use $a$ modulo fourth powers.

For a candidate $\mathbb F_{p^d}$, put $M=p^n-1$, $m=p^d-1$ and $h=\gcd(\text{twist order},M/m)$. The image subgroup $\mathbb F_{p^d}^{*}(K^{*})^{\text{twist order}}$ is generated by the $h$-th power of a cyclic generator. Thus the twist parameter lies in it exactly when its $M/h$ power is one. Each row retains all candidate power witnesses, including rejected smaller fields.

Four $p=7$ records have $j$ degree three. Three have a $K$-isomorphic model over $\mathbb F_{7^3}$, and one requires $\mathbb F_{7^6}$. Only one already displays coefficients over $\mathbb F_{7^3}$:
\begin{center}\begin{tabular}{lrrr}\toprule Exact model slug & $j$ & Coefficients & Isomorphic model\\\midrule
\texttt{icv1-fp3k6-tm110-5f712c89} &3&6&3\\
\texttt{icv1-fp3k6-t202-39e02f24} &3&6&6\\\bottomrule\end{tabular}\end{center}
These are model-definition properties. Selected-point and subgroup containment remain unspecified.
\section{Certificates and verification}
Thirteen portable identity certificates cover the cubic discriminant, cleared normalized-$j$ and twist coefficients, translated Frobenius norm and eight unit-group index degrees. Each passes 64 saved and 128 fresh evaluations; native expansion additionally checks exact integer coefficients. IDs and expression trees are in \texttt{identity\_certificates.json}; the formal-proof pointer remains null.

Native arithmetic and PARI agree on all 814 property rows. Audit run 2 completes under its 300-second cap with empty stderr. The initial stack-setting initialization failure, its zero-row receipt and executed sources are retained. Point counts are reused from the independently verified parent records. Ten $p=7$ records have 13-smooth ambient order; residual cofactors are unclassified, and selected-subgroup costs are unmeasured.

The property browser is a sidecar. The earlier catalogue study and two-branch evidence remain sealed, and numerical IC ratio/census graphs retain their values: this round adds structural classification rather than a new paired timing measurement. Reproduction, source hashes, identity IDs and raw arithmetic evidence are linked in the Markdown report and validation record.
\section{Sources}
\href{https://math.mit.edu/classes/18.783/2023/LectureNotes13.pdf}{Sutherland: ordinary/supersingular curves and endomorphism orders}.\par
\href{https://math.mit.edu/classes/18.783/2023/LectureNotes22.pdf}{Sutherland: isogeny volcanoes}.\par
\href{https://pari.math.u-bordeaux.fr/dochtml/html-stable/Arithmetic_functions.html}{PARI: arithmetic and fundamental discriminants}.\par
\href{https://pari.math.u-bordeaux.fr/dochtml/html/Polynomials_and_power_series.html}{PARI: finite-field polynomial roots}.\par
\href{https://github.com/aburan28/crypto/blob/862f35365af35cdbc29a2cd9b6ab63037afb668c/research/iso1_weak_classes_20261007/two_branch_20261009/REPORT.md}{Parent two-branch theorem and frozen model evidence}.
\end{document}
"#;
pub fn check(repo: &Path, data: &Value) {
    assert_eq!(
        fs::read_to_string(repo.join("docs/curves/weak-families/properties.html")).unwrap(),
        html(data)
    );
}
const HTML: &str = r##"<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Certified curve properties</title>
<style>:root{color-scheme:light dark;--bg:#f6f8fb;--ink:#17304c;--card:white;--line:#dce3ed;--muted:#52667f;--accent:#245aaa}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 system-ui,sans-serif}main{max-width:1300px;margin:auto;padding:28px 24px}h1{font-size:34px;margin:0 0 10px}h2{font-size:22px}a{color:var(--accent)}.intro{max-width:1000px}.small{font-size:13px;color:var(--muted)}.controls{display:flex;gap:9px;flex-wrap:wrap;margin:20px 0}input,select,button{font:inherit;border:1px solid var(--line);border-radius:7px;background:var(--card);color:var(--ink);padding:9px 11px}input{flex:1;min-width:220px}.scroll{overflow:auto;border:1px solid var(--line);border-radius:9px;background:var(--card)}table{width:100%;min-width:1000px;border-collapse:collapse}td,th{padding:11px 12px;text-align:left;vertical-align:top;border-bottom:1px solid var(--line)}.note,.detail{background:var(--card);border:1px solid var(--line);padding:15px 18px;border-radius:8px;overflow-wrap:anywhere}.detail{margin:20px 0}.hide{display:none}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}code{font:12px/1.5 ui-monospace,monospace;overflow-wrap:anywhere}.box{margin:12px 0}.badge{font-size:12px;background:var(--bg);padding:3px 7px;border-radius:4px}.controls>*{max-width:100%}@media(max-width:700px){main{padding:20px 13px}h1{font-size:27px}input{min-width:100%;width:100%}}@media(prefers-color-scheme:dark){:root{--bg:#111b29;--ink:#e4edf8;--card:#182638;--line:#34455c;--muted:#b1bfd3;--accent:#8cb9ff}}</style>
<main><h1>Certified curve properties</h1><p class="intro">Frobenius conductor depth, rational torsion and three fields of definition, checked in native arithmetic and PARI. Family membership and class support retain their existing evidence labels.</p><p class="small" id="stats"></p><p><a href="index.html">Family catalogue and original class table</a> · <a href="properties.json">Property records and witnesses</a> · <a href="../../../research/weak_curve_properties_20261010/REPORT.md">Methods and evidence</a></p>
<div class="note">The degree of j, the degree of the displayed coefficients and the degree of an isomorphic model answer separate questions. A smaller model field leaves selected-point and subgroup containment unspecified. The conductor column is for the Frobenius order; the full endomorphism conductor stays unknown. Local 2-levels are exact where torsion certifies them, with bounds elsewhere.</div>
<div class="controls"><input id="query" type="search" placeholder="Search a model ID or trace" aria-label="Search model ID or trace"><select id="cohort" aria-label="Inventory cohort"><option value="">All inventories</option><option value="p7_replay_control">p7 replay and controls</option><option value="higher_degree_control">Higher-degree controls</option><option value="independent_large">Independent large sources</option></select><select id="property" aria-label="Property filter"><option value="">All properties</option><option value="full2">Full rational 2-torsion</option><option value="full4">Full rational 4-torsion</option><option value="jproper">j in a proper subfield</option><option value="literalproper">Displayed coefficients in a proper subfield</option><option value="isoproper">Isomorphic model in a proper subfield</option><option value="jblocked">j descends; isomorphic model stays in full field</option><option value="smooth">Ambient order is 13-smooth</option><option value="ss">Supersingular</option><option value="exactlevel">Certified ordinary 2-level</option><option value="floor">Certified ordinary 2-floor</option></select><select id="depth" aria-label="Frobenius conductor depth"><option value="">All Frobenius 2-depths</option><option value="ge2">Depth at least 2</option></select><select id="splitting" aria-label="Splitting of two"><option value="">All 2-splitting types</option><option value="split">Split</option><option value="inert">Inert</option><option value="ramified">Ramified</option></select></div>
<p id="count" class="small"></p><section id="detail" class="detail hide"></section><div class="scroll"><table><thead><tr><th>Exact model</th><th>Field bits / degree</th><th>Frobenius 2-depth / splitting</th><th>Endomorphism 2-level</th><th>Rational nonzero 2-points</th><th>Full 4-torsion</th><th>Degrees: j / coefficients / isomorphic model</th></tr></thead><tbody id="rows"></tbody></table></div><div class="controls"><button id="previous">Previous</button><span id="page"></span><button id="next">Next</button></div><noscript><p>Canonical property records and the report are available in the links above. Filtering requires JavaScript.</p></noscript></main>
<script type="application/json" id="data">__DATA__</script><script>
const data=JSON.parse(document.getElementById('data').textContent),$=id=>document.getElementById(id);let page=0,filtered=[];const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));const pretty=s=>String(s??'').replaceAll('_',' ');$('stats').textContent=`${data.records.length} exact model records · native and independent PARI verification · unchanged parent family/class labels`;
const depths=[...new Set(data.records.map(r=>r.properties.cm.frobenius_conductor_v2).filter(d=>d!==null))].sort((a,b)=>a-b);$('depth').innerHTML+=depths.map(d=>`<option value="${d}">Depth ${d}</option>`).join('');
function draw(){const query=$('query').value.toLowerCase(),cohort=$('cohort').value,property=$('property').value,depth=$('depth').value,splitting=$('splitting').value;filtered=data.records.filter(r=>{const p=r.properties,d=p.cm.frobenius_conductor_v2;return(!query||(r.icv1+' '+r.trace+' '+r.slug).toLowerCase().includes(query))&&(!cohort||r.cohort===cohort)&&(!splitting||p.cm.two_splitting===splitting)&&(!depth||(depth==='ge2'?d!==null&&d>=2:String(d)===depth))&&(!property||({full2:p.rational_2_nonidentity_points===3,full4:p.full_rational_4_torsion,jproper:p.j_in_proper_subfield,literalproper:p.displayed_model_in_proper_subfield,isoproper:p.isomorphic_model_in_proper_subfield,jblocked:p.j_in_proper_subfield&&!p.isomorphic_model_in_proper_subfield,smooth:p.ambient_order_13_smooth,ss:p.supersingular,exactlevel:p.cm.endomorphism_conductor_v2!==null,floor:p.cm.endomorphism_conductor_v2!==null&&p.cm.endomorphism_conductor_v2===p.cm.frobenius_conductor_v2}[property]));});page=Math.max(0,Math.min(page,Math.max(0,Math.ceil(filtered.length/50)-1)));$('count').textContent=`${filtered.length} model records; counts describe the selected inventory`;$('rows').innerHTML=filtered.slice(page*50,(page+1)*50).map(r=>{const p=r.properties;return `<tr><td><button class="model" data-id="${esc(r.icv1)}">${esc(r.slug)}</button><br><span class="small">${esc(pretty(r.cohort))}</span></td><td>${r.field.bit_length} / ${r.field.degree}</td><td>${p.cm.frobenius_conductor_v2??'Not applicable'} / ${esc(pretty(p.cm.two_splitting))}</td><td>${p.cm.endomorphism_conductor_v2??(p.cm.endomorphism_conductor_v2_bounds?`0–${p.cm.endomorphism_conductor_v2_bounds[1]} (bounded)`:"Not applicable")}</td><td>${p.rational_2_nonidentity_points}</td><td>${p.full_rational_4_torsion?'Yes':'No'}</td><td>${p.j_minimum_degree} / ${p.coefficient_minimum_degree} / ${p.isomorphic_model_minimum_degree}</td></tr>`;}).join('');$('page').textContent=`Page ${page+1} / ${Math.max(1,Math.ceil(filtered.length/50))}`;$('previous').disabled=page===0;$('next').disabled=(page+1)*50>=filtered.length;document.querySelectorAll('.model').forEach(b=>b.onclick=()=>detail(data.records.find(r=>r.icv1===b.dataset.id)));}
function detail(r){const p=r.properties;$('detail').classList.remove('hide');$('detail').innerHTML=`<h2>${esc(r.slug)}</h2><p><code>${esc(r.icv1)}</code></p><p><strong>Trace:</strong> ${esc(r.trace)}<br><strong>Order:</strong> ${esc(r.order)}<br><strong>Field degrees:</strong> j ${p.j_minimum_degree}, displayed coefficients ${p.coefficient_minimum_degree}, isomorphic model ${p.isomorphic_model_minimum_degree}<br><strong>Selected points and subgroup:</strong> containment unspecified; comparative cost unmeasured.</p><p><strong>Frobenius conductor:</strong> 2-depth ${p.cm.frobenius_conductor_v2??'not applicable'}; <strong>endomorphism 2-depth:</strong> ${p.cm.endomorphism_conductor_v2??"not determined"}; full conductor unknown${p.cm.endomorphism_conductor_v2_bounds?`, permitted 2-depth from ${p.cm.endomorphism_conductor_v2_bounds[0]} to ${p.cm.endomorphism_conductor_v2_bounds[1]}`:''}.</p><details><summary>Existing family and class evidence</summary><pre>${esc(JSON.stringify(r.family_labels,null,2))}</pre></details><details><summary>Property-derived family checks</summary><pre>${esc(JSON.stringify(r.family_property_checks,null,2))}</pre></details><details><summary>Exact properties and arithmetic witnesses</summary><pre>${esc(JSON.stringify({properties:p,witnesses:r.witnesses},null,2))}</pre></details>`;$('detail').scrollIntoView({block:'nearest'});}
for(const id of ['query','cohort','property','depth','splitting'])$(id).addEventListener(id==='query'?'input':'change',()=>{page=0;draw()});$('previous').onclick=()=>{page--;draw()};$('next').onclick=()=>{page++;draw()};draw();
</script></html>"##;
