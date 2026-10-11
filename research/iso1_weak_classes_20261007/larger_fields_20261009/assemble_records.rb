# Assemble documentary identities from frozen GP vectors; perform no curve arithmetic.
require 'json'
require 'digest'
require 'csv'
root=__dir__
receipts=CSV.read(File.join(root,'evidence_run1/receipts.tsv'),headers:true,col_sep:"\t")
sizes=File.readlines(File.join(root,'field_sizes.txt')).to_h do |line|
  p,n,q,bits,c,logc=line.strip.split('|'); [[p,n],{q:q,log2_q:bits,orbit_calls:c,log10_calls:logc}]
end
raise 'wrong cell count' unless receipts.length==34 && receipts.all?{|r|r['status']=='COMPLETE'}
def identity(p,k,mod,a,b,j,order,trace)
  field="fpk-#{p}-#{k}-#{Digest::SHA256.hexdigest("fpk-modulus:#{p}:#{mod.join(',')}")[0,8]}"
  model={'a'=>a.map(&:to_s),'b'=>b.map(&:to_s),'field'=>field,'form'=>'y^2=x^3+a*x+b','k'=>k.to_s,'modulus'=>mod.map(&:to_s),'p'=>p.to_s,'v'=>'1'}
  canonical=JSON.generate(model.sort.to_h)
  h=Digest::SHA256.hexdigest(canonical)
  signed=trace<0 ? "tm#{-trace}" : "t#{trace}"
  {'icv1'=>"ICV1:#{field}:#{trace}:#{order}:#{j.join(',')}:unk:unk:r:#{h[0,12]}",
   'slug'=>"icv1-fp#{p.bit_length}k#{k}-#{signed}-#{h[0,8]}",'model_json'=>canonical,'model_sha256'=>h,
   'trace'=>trace.to_s,'group_order'=>order.to_s,'j_coefficients'=>j.map(&:to_s),
   'endomorphism_order_conductor'=>nil,'frobenius_order_conductor'=>nil,'volcano_level'=>nil,'subgroup_order'=>nil}
end
pairs=[]
receipts.select{|r|r['mode']=='count'}.each do |r|
  p=r['p'].to_i;n=r['odd_degree'].to_i;k=2*n
  raw=File.read(File.join(root,"evidence_run1/p#{p}_n#{n}_count.stdout"))
  field=raw.lines.find{|l|l.start_with?('FIELD|')}.strip.split('|')
  mod=JSON.parse(field[5]);q=Integer(field[4])
  g=raw.lines.find{|l|l.start_with?('GEOMETRY|')}.strip.split('|')
  geometry=File.read(File.join(root,"evidence_run1/p#{p}_n#{n}_geometry.stdout"))
  geometry_field=geometry.lines.find{|l|l.start_with?('FIELD|')}.strip.split('|')
  geometry_first=geometry.lines.find{|l|l.start_with?("GEOMETRY|#{p}|1|")}.strip.split('|')
  raise 'geometry/count field mismatch' unless field == geometry_field
  raise 'geometry/count lambda mismatch' unless JSON.parse(g[4]) == JSON.parse(geometry_first[4])
  row=raw.lines.find{|l|l.start_with?('COUNT|')}.strip.split('|')
  order=Integer(row[5]);trace=Integer(row[6]);raise 'bad count metadata' unless q+1-order==trace && order%16==0 && trace*trace<=4*q
  vectors=row[10..15].map{|s|JSON.parse(s)}
  raise 'noncanonical coefficient vector' unless [mod,*vectors].all?{|v|v.length==k&&v.all?{|x|x>=0&&x<p}}
  source=identity(p,k,mod,*vectors[0..2],order,trace)
  target=identity(p,k,mod,*vectors[3..5],order,trace)
  size=sizes.fetch([p.to_s,n.to_s])
  pairs << {'characteristic'=>p,'base_degree'=>2,'odd_degree'=>n,'field_degree'=>k,
    'field_cardinality'=>q.to_s,'log2_field_cardinality'=>size[:log2_q],
    'absolute_orbit_calls'=>size[:orbit_calls],'log10_absolute_orbit_calls'=>size[:log10_calls],
    'modulus_low_coefficients'=>mod,'lambda_coefficients'=>JSON.parse(g[4]),
    'ordinary'=>row[9]=='1','source_count_ms'=>row[7].to_i,'target_count_ms'=>row[8].to_i,
    'source'=>source,'target'=>target,'isogeny'=>{'degree'=>2,'source'=>source['icv1'],'target'=>target['icv1'],
      'kernel_original_model'=>[[0,0]],'original_map'=>'(x,y) -> (x+a+b/x, y*(1-b/x^2)); a=-(1+lambda), b=lambda',
      'short_model_coordinate_changes'=>'X=x+a2/3 on each original model; Y=y',
      'status'=>'explicit_point_replay_and_equal_independent_cardinalities','identity_record'=>'IDC1he02ab9a5cf019a56'},
    'verification'=>{'norm_one'=>true,'fourth_root'=>true,'two_independent_target_halves'=>true,'full_rational_target_4_torsion'=>true,'point_count_pair_equal'=>true,'hasse'=>true,'cardinality_mod16'=>0}}
end
summary={'schema'=>'iso1-larger-field-controls/v1','date'=>'2026-10-09','geometry_controls'=>49,'paired_counts'=>pairs.length,'field_cells'=>17,'input_law'=>'uniform nonidentity norm-one lambda from random(z)^(p^2-1), per-curve seed 20261009+i','curve_identity_spec'=>'docs/curves/ICV1.md','prime_degree_orbit_identity'=>'IDC1h90d58cc0e0c48fe3','records'=>pairs}
File.write(File.join(root,'curve_records.json'),JSON.pretty_generate(summary)+"\n")
headers=%w[p odd_degree field_degree log2_field_size source_order source_trace source_count_ms target_count_ms orbit_calls source_slug target_slug]
CSV.open(File.join(root,'summary.csv'),'w'){|csv|csv<<headers;pairs.each{|r|csv<<[r['characteristic'],r['odd_degree'],r['field_degree'],r['log2_field_cardinality'],r['source']['group_order'],r['source']['trace'],r['source_count_ms'],r['target_count_ms'],r['absolute_orbit_calls'],r['source']['slug'],r['target']['slug']]}}
puts "#{pairs.length} paired counts; 49 geometry controls; #{pairs.count{|r|r['ordinary']}} ordinary; #{pairs.length*2} ICV1 model identities"
