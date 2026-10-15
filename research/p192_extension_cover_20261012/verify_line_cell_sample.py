"""Custody and selection-law audit for the exploratory P-192 cell sample."""

import hashlib
import itertools
import json
import random
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
v=json.loads((HERE/"type_i_line_cell_enrichment_sample.json").read_text())
assert v["status"]=="bounded_p192_line_support_vs_uniform_character_cell_sample"
assert v["source_sha256"]==sha(HERE/"line_cell_enrichment_sample.py")
assert v["base_receipt_sha256"]==sha(HERE/"type_i_line_guided_collector.json")
assert v["support_receipt_sha256"]==sha(HERE/"type_i_line_three_sum_coverage.json")
assert v["runtime_sha256"]==sha(HERE/"type_i_line_cell_sample_runtime_info.json")
assert v["support_population"]==477191
assert v["support_sample_without_replacement"]==v["uniform_scalar_sample_with_replacement"]==512
assert len(v["support_rows"])==len(v["uniform_rows"])==512
rng=random.Random(v["seed"])
indices=sorted(rng.sample(range(477191),512))
assert [row["triple_index"] for row in v["support_rows"]]==indices
wanted=set(indices)
selected=[(index,[i,j,k]) for index,(i,j,k) in
          enumerate(itertools.combinations_with_replacement(range(141),3))
          if index in wanted]
assert [(row["triple_index"],row["indices"]) for row in v["support_rows"]]==selected
r=int(v["subgroup_order"])
for index,row in enumerate(v["uniform_rows"]):
    assert row["trial_index"]==index
    assert int(row["scalar"])==rng.randrange(1,r)
valid={"++","+-","-+","--","exceptional"}
for label,rows in (("support",v["support_rows"]),("uniform",v["uniform_rows"])):
    assert all(row["cell"] in valid for row in rows)
    assert Counter(row["cell"] for row in rows)=={
        k:int(val) for k,val in v[label+"_cells"].items() if val}
assert v["support_cells"]=={"++":131,"+-":103,"-+":140,"--":138,"exceptional":0}
assert v["uniform_cells"]=={"++":125,"+-":130,"-+":133,"--":124,"exceptional":0}
print("LINE_CELL_SAMPLE_RECEIPT_PASS",json.dumps({
    "support_cells":v["support_cells"],"uniform_cells":v["uniform_cells"]},
    sort_keys=True))
