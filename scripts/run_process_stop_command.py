"""Hard60s supervision and immutable command receipts for this work package."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--label',required=True)
    ap.add_argument('command',nargs=argparse.REMAINDER)
    a=ap.parse_args()
    if not a.label.replace('_','').isalnum():raise ValueError('simple fresh label required')
    cmd=a.command[1:] if a.command[:1]==['--'] else a.command
    if not cmd or cmd[0] not in ('scripts/check_process_stop.py','scripts/probe_hht_stop.py'):
        raise ValueError('only registered component gate or local scout command')
    directory=ROOT/'logs/process_stop_20261003';directory.mkdir(exist_ok=True)
    receipt=directory/(a.label+'_supervision.json');log=directory/(a.label+'.log')
    if receipt.exists() or log.exists():raise FileExistsError('never overwrite command evidence')
    command=[sys.executable,'-B']+cmd
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    started=time.perf_counter()
    with log.open('x',encoding='utf-8') as stream:
        process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
        try:code=process.wait(timeout=60.)
        except subprocess.TimeoutExpired:
            process.kill();process.wait();stream.write('\nEXTERNAL60S_TIMEOUT\n');code=124
    elapsed=time.perf_counter()-started
    data=dict(command=command,return_code=code,wall_seconds=elapsed,external_limit_seconds=60.,
        log=log.name,log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
        supervisor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    receipt.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())
