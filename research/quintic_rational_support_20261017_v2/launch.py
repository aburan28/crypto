"""Existing user-required legacy scheduler adapter; all scientific work is native."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('quintic-rational-support'))
os.execv(p,[p,*sys.argv[1:]])
