"""Use .envs/pysph/python.exe; all outputs are explicitly DEBUG."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.application import main

if __name__ == '__main__':
    raise SystemExit(main())
