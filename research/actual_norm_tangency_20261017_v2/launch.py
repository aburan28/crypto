"""Compatibility with the user-required scheduler; native science only."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('actual-norm-tangency'))
os.execv(p,[p,*sys.argv[1:]])
