"""Compatibility with the user-required scheduler; compiler and tests native."""
import os,sys
from pathlib import Path
p=str(Path(__file__).with_name('direct-support-karatsuba'))
os.execv(p,[p,*sys.argv[1:]])
