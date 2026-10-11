# Documentary context only: preserve every quantitative progress point.
require 'json'
repo=ARGV.fetch(0)
path=File.join(repo,'docs/ic/progress-timeline.json')
raw=File.read(path)
data=JSON.parse(raw)
cover=data.fetch('series').find{|s|s.fetch('id')=='cover-fp6'}
raise 'missing cover context' unless cover
addition=' The independently sampled 192-252-bit population on 2026-10-09 point-counts 512 full-2-torsion curves in ten fields of degrees 6, 10 and 14. The conductor condition admits 323/512 (63.0859%, pooled Wilson 95% 58.8229-67.1540%); all 323 have a verified full-4 representative within one rational 2-isogeny after choosing the trace sign. A new theorem proves this exact torsion equivalence. Bounded degree-2/3 searches test 32,484 vertices and leave all 323 admitted weak-class labels unresolved (219 restricted component closures, 104 vertex caps); all ten separately constructed positive controls recover a weak endpoint. The direct weak-model density is bounded by 3/(p^2-1). A recorded exact CM preflight obtains f_pi=8 on a 192-bit source, then polclass fails with an integer-conversion overflow. These curve-weighted admissions do not replace the small-prime class precision or alter IC/rho ratio points.'
unless cover.fetch('note').include?('independently sampled 192-252-bit population')
  old=JSON.generate(cover.fetch('note'))
  replacement=JSON.generate(cover.fetch('note')+addition)
  raise 'ambiguous source note' unless raw.scan(old).length==1
  raw=raw.sub(old,replacement)
  File.write(path,raw)
end
correction=' Scope correction, 2026-10-09: the census and 63.1% admission concern the cubic norm-one h=x branch. Joux-Vitse also permit nonsplit quadratic h. Three exhaustive quadratic controls verify traces -38, -10 and 610 over F_(7^6), with orders 117688, 117660 and 117040. Hence the cubic zero at +/-610 and depth-one exclusions are not full-family exclusions; the next classification must cover both branches.'
['cover-fp6','cover-walk-stage'].each do |id|
  current=JSON.parse(raw).fetch('series').find{|v|v.fetch('id')==id}
  unless current.fetch('note').include?('Scope correction, 2026-10-09')
    old=JSON.generate(current.fetch('note'))
    raw=raw.sub(old,JSON.generate(current.fetch('note')+correction))
  end
end
File.write(path,raw)

html_path=File.join(repo,'docs/index-calculus-scoreboard.html')
html=File.read(html_path)
pattern=/(<script type="application\/json" id="progress-data">\n).*?(\n<\/script>)/m
raise 'progress-data block absent' unless html.match?(pattern)
html=html.sub(pattern){$1+raw.strip+$2}
unless html.include?('<strong>192-252-bit independent population')
  marker='<p class="dash-section-intro"><strong>Larger norm-one fields, 9 October 2026:</strong>'
  raise 'ISO-1 context marker absent' unless html.include?(marker)
  paragraph='<p class="dash-section-intro"><strong>192-252-bit independent population, 9 October 2026:</strong> 323 of 512 independently sampled full-2-torsion curves pass the conductor condition (63.1%; pooled Wilson 95% 58.8-67.2%). All admitted curves have a verified full-4 representative within one rational 2-isogeny; the new theorem proves the equivalence. Degree-2/3 searches test 32,484 vertices and leave all 323 admitted weak-class labels unresolved; the ten constructed positive controls all succeed. The report retains exact fields, moduli, model identities, caps, raw receipts, and the CM preflight overflow. <a href="https://github.com/aburan28/crypto/blob/main/research/iso1_weak_classes_20261007/large_population_20261009/REPORT.md">Population measurement and proofs</a>.</p>'
  html=html.sub(marker,paragraph+"\n"+marker)
end
unless html.include?('<strong>Cubic and quadratic branch scope:')
  marker='<p class="dash-section-intro"><strong>192-252-bit independent population'
  paragraph='<p class="dash-section-intro"><strong>Cubic and quadratic branch scope:</strong> the conductor theorem, census zeros and 63.1% admission describe the norm-one cubic branch. The Joux-Vitse paper also includes nonsplit quadratic h. An exact quadratic model at p = 7 has trace 610 and 117,040 points, although that class is empty in the cubic census; two other quartics have conductor depth one. <a href="https://github.com/aburan28/crypto/blob/main/research/iso1_weak_classes_20261007/large_population_20261009/PRIOR_WORK.md">Attribution and independently verified scope correction</a>.</p>'
  html=html.sub(marker,paragraph+"\n"+marker)
end
html=html.sub('has conductor 72 and no weak representative','has conductor 72 and no cubic norm-one representative')
File.write(html_path,html)

embedded=JSON.parse(html.match(pattern)[0].sub(/\A.*?\n/m,'').sub(/\n<\/script>\z/,''))
raise 'source/embedded drift' unless embedded==JSON.parse(raw)
puts 'canonical context and embedded copy agree; ratio arrays unchanged'
