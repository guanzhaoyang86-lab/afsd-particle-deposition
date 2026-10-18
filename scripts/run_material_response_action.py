"""Bound one saved-data action process; no solver calls by this supervisor."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--label',required=True);ap.add_argument('command',nargs=argparse.REMAINDER)
args=ap.parse_args();command=args.command
if command and command[0]=='--':command=command[1:]
if not command or command[0]!='scripts/probe_material_response_action.py':raise ValueError('registered action only')
if not args.label.replace('_','').isalnum():raise ValueError('simple label')
out=ROOT/'logs/process_material_response_20261004';out.mkdir(exist_ok=True)
log=out/(args.label+'.log');receipt=out/(args.label+'_supervision.json')
if log.exists() or receipt.exists():raise FileExistsError('fresh command label')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
cmd=[sys.executable,'-B']+command;start=time.perf_counter()
with log.open('x',encoding='utf-8') as stream:
    proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
    try:code=proc.wait(timeout=60.)
    except subprocess.TimeoutExpired:
        proc.kill();proc.wait();stream.write('\nEXTERNAL60S_TIMEOUT\n');code=124
r=dict(command=cmd,return_code=code,wall_seconds=time.perf_counter()-start,external_limit_seconds=60,
       new_physical_evaluations=0,log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
       supervisor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
with receipt.open('x',encoding='utf-8') as f:json.dump(r,f,indent=2,allow_nan=False)
print(json.dumps(r),flush=True)
raise SystemExit(code)

