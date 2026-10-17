import os,sys
from pathlib import Path
p=str(Path(__file__).with_name("hybrid-norm-search"))
os.execv(p,[p,*sys.argv[1:]])
