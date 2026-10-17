"""Compatibility with user-required scheduler; compiler and oracles native."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('crt-support-circuit'))
os.execv(p,[p,*sys.argv[1:]])
