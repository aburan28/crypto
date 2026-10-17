"""Compatibility with the user-required scheduler; all science native."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('odd-interval-roots'))
os.execv(p,[p,*sys.argv[1:]])
