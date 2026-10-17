"""User-required existing scheduler adapter; native field/SAT work only."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('full-binary-norm-sat'))
os.execv(p,[p,*sys.argv[1:]])
