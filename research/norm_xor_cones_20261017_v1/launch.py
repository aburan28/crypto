"""Existing user-required scheduler adapter; all mathematical work native."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('norm-xor-cones'))
os.execv(p,[p,*sys.argv[1:]])
