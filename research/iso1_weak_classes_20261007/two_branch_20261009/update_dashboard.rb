# Documentary context; retain all existing quantitative IC/rho points.
require 'json'
repo=ARGV.fetch(0)
path=File.join(repo,'docs/ic/progress-timeline.json')
raw=File.read(path)
original=JSON.parse(raw)
note=' The two-branch supplement of 2026-10-09 proves complete cubic and nonsplit-quadratic torus normalizations and an exact CM class criterion with distinct conductor selections. At p37 their union covers 48630/49284 ordinary classes (cubic 24352, quadratic 48164, overlap 23886, exact zeros 654). Independent all-x replay verifies 409 p7 models in 48118441 evaluations; all 63 positive fixtures among 64 independent p7 sources have explicit verified routes after broader-prime follow-ups. The expanded 512-source 192-252-bit search visits 46356 vertices and evaluates 129600 degree-2/3 edges; its 383 restricted closures and 129 vertex caps leave all large-field class labels unresolved. Quadratic-base Hilbert 90 validates 16 degree-10/14 controls, eight with noncoprime section gcd5; the specified compositum cover genera are 25 and 161. Finite-field CM decision degree is at most18 per factor, with CM support generation still a separate cost.'
['cover-fp6','cover-walk-stage'].each do |id|
  row=JSON.parse(raw).fetch('series').find{|s|s.fetch('id')==id}
  raise 'missing series' unless row
  unless row.fetch('note').include?('The two-branch supplement of 2026-10-09')
    old=JSON.generate(row.fetch('note'))
    raise 'ambiguous source note' unless raw.scan(old).length==1
    raw=raw.sub(old,JSON.generate(row.fetch('note')+note))
  end
end
updated=JSON.parse(raw)
original.fetch('series').zip(updated.fetch('series')).each do |a,b|
  raise 'quantitative series changed' unless a.reject{|k,_|k=='note'}==b.reject{|k,_|k=='note'}
end
File.write(path,raw)
html_path=File.join(repo,'docs/index-calculus-scoreboard.html')
html=File.read(html_path)
pattern=/(<script type="application\/json" id="progress-data">\n).*?(\n<\/script>)/m
raise 'embedded source missing' unless html.match?(pattern)
html=html.sub(pattern){$1+raw.strip+$2}
unless html.include?('<strong>Exact support of both genus-3 branches:')
  marker='<p class="dash-section-intro"><strong>Cubic and quadratic branch scope:'
  raise 'scope context missing' unless html.include?(marker)
  paragraph='<p class="dash-section-intro"><strong>Exact support of both genus-3 branches:</strong> the combined torus/CM criterion covers 48,630 of 49,284 ordinary trace classes at p = 37; 654 are exact zeros for this hyperelliptic family. All 63 positive fixtures among 64 independent p = 7 sources have verified routes. The expanded 512-source 192-252-bit search visits 46,356 vertices and leaves its large-field class labels unresolved within the recorded caps. The degree-10/14 quadratic reconstruction passes 16 controls; their specified compositum covers have genera 25 and 161. <a href="https://github.com/aburan28/crypto/blob/main/research/iso1_weak_classes_20261007/two_branch_20261009/REPORT.md">Exact decisions, construction receipts and full proofs</a>.</p>'
  html=html.sub(marker,paragraph+"\n"+marker)
end
File.write(html_path,html)
embedded=JSON.parse(html.match(pattern)[0].sub(/\A.*?\n/m,'').sub(/\n<\/script>\z/,''))
raise 'embedded JSON drift' unless embedded==updated
puts 'canonical context and rendered copy agree; all quantitative IC/rho arrays preserved'
