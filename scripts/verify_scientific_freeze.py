"""Read-only evidence/physics freeze verification; never updates the manifest."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def verify(path):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    errors = []
    for relative, expected in data['files_sha256'].items():
        candidate = (ROOT / relative).resolve()
        if not candidate.is_relative_to(ROOT) or not candidate.is_file():
            errors.append(f'missing/unsafe path: {relative}')
            continue
        actual = hashlib.sha256(candidate.read_bytes()).hexdigest()
        if actual != expected:
            errors.append(f'SHA256 MISMATCH: {relative}')
    pin = data['backend_commit']
    repository = ROOT/'third_party/PySPH'
    actual_pin = subprocess.check_output(['git', '-C', str(repository), 'rev-parse', 'HEAD'], text=True).strip()
    dirty = subprocess.check_output(['git', '-C', str(repository), 'status', '--porcelain', '--untracked-files=no'], text=True).strip()
    if actual_pin != pin or dirty:
        errors.append('upstream pin or tracked-source cleanliness mismatch')
    print(json.dumps({'status': 'FAILED' if errors else 'PASS',
                      'files_checked': len(data['files_sha256']),
                      'errors': errors, 'production_authorized': False,
                      'scope': 'file integrity only; NOT scientific acceptance'}, indent=2))
    return 1 if errors else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', default=str(ROOT/'docs/SCIENTIFIC_FREEZE.json'))
    raise SystemExit(verify(parser.parse_args().manifest))
