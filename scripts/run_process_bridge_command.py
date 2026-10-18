"""Fresh evidence: preflight60s, exactly-four-map bridge120s external limits."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--label',required=True)
    ap.add_argument('command',nargs=argparse.REMAINDER);args=ap.parse_args()
    if not args.label.replace('_','').isalnum():raise ValueError('simple label required')
    cmd=args.command[1:] if args.command[:1]==['--'] else args.command
    if not cmd or cmd[0] not in ('scripts/check_process_bridge.py','scripts/probe_afsd_process_bridge.py'):
        raise ValueError('only bridge preflight or registered measurement')
    folder=ROOT/'logs/process_bridge_20261003';folder.mkdir(parents=True,exist_ok=True)
    log=folder/(args.label+'.log');receipt=folder/(args.label+'_supervision.json')
    if log.exists() or receipt.exists():raise FileExistsError('never overwrite evidence')
    command=[sys.executable,'-B']+cmd;started=time.perf_counter()
    limit=120. if cmd[0]=='scripts/probe_afsd_process_bridge.py' else 60.
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    with log.open('x',encoding='utf-8') as stream:
        proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
        try:code=proc.wait(timeout=limit)
        except subprocess.TimeoutExpired:
            if os.name=='nt':
                subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=5,check=False)
            if proc.poll() is None:proc.kill()
            proc.wait();stream.write('\nEXTERNAL_WALL_TIMEOUT\n');code=124
    data=dict(command=command,return_code=code,wall_seconds=time.perf_counter()-started,
        external_limit_seconds=limit,new_published_steps=0,production_authorized=False,
        log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
        supervisor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with receipt.open('x',encoding='utf-8') as stream:json.dump(data,stream,indent=2,allow_nan=False)
    print(json.dumps(data),flush=True);return code


if __name__=='__main__':raise SystemExit(main())
