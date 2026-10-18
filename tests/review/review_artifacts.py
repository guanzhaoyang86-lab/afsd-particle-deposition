"""Keep newly generated benchmark reports out of frozen historical evidence."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def artifact_path(name):
    directory=Path(os.environ.get('AFSD_REVIEW_ARTIFACT_DIR',str(ROOT/'outputs/review_rechecks')))
    directory.mkdir(parents=True,exist_ok=True)
    return directory/name
