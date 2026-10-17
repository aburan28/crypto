"""Existing user-required legacy scheduler adapter; scientific logic is native."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('quintic-fixed-u-cells'))
os.execv(p,[p,*sys.argv[1:]])
