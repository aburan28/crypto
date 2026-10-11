from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
