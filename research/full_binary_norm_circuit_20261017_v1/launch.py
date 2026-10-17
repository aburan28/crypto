"""Compatibility with user-required scheduler; circuit and oracles native."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('full-binary-norm-circuit'))
os.execv(p,[p,*sys.argv[1:]])
