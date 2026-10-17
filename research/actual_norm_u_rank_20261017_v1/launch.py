"""Compatibility with the user-required scheduler; native scientific work only."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('actual-norm-u-rank'))
os.execv(p,[p,*sys.argv[1:]])
