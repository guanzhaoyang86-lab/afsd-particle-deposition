"""Two-step continuation of a real release checkpoint; not a new physical run."""
from pathlib import Path
import argparse
import json
import time
from probe_axial_feed import Block
from probe_rotating_deformable_block import jsonable
from inspect_rotating_deformable_block import load,file_hash,require,close


def check(directory):
    start=time.perf_counter();directory=Path(directory);m=load(directory/'manifest.json')
    require(m['accepted_time_s']==8e-7 and m['production_authorized'] is False,'complete registered source case')
    h=load(directory/'history.json');path=directory/'release_accepted.npz'
    require(file_hash(path)==m['artifact_sha256'][path.name],'checkpoint hash')
    require(file_hash(directory/'history.json')==m['artifact_sha256']['history.json'],'history hash')
    b=Block.from_checkpoint(path);require(b.step in (80,160),'release step')
    close(b.time,4e-7,1e-20,'release time');first=b.step;tested=[];differences={}
    for offset in range(3):
        r=jsonable(b.record());ref=h[b.step]
        for key in ('particle_id','tracer_type','as_rs_tag'):
            require(r[key]==ref[key],'restart identity '+key)
        for field,tol in [('position_m',1e-12),('velocity_m_s',1e-12),('density_kg_m3',1e-12),
            ('deviatoric_stress_Pa',1e-6),('ep',1e-12)]:
            diff=close(r[field],ref[field],tol,'restarted '+field)
            differences[field]=max(differences.get(field,0.),diff)
        for field,value in ref['cumulative'].items():close(r['cumulative'][field],value,1e-15,'restarted ledger '+field)
        tested.append(b.step)
        if offset<2:b.advance()
    return dict(status='PASS_NONZERO_RELEASE_CHECKPOINT_TWO_STEPS_NOT_AFSD',production_authorized=False,
        manifest_sha256=file_hash(directory/'manifest.json'),checkpoint_sha256=file_hash(path),
        script_sha256=file_hash(Path(__file__)),source_sha256=m['source_sha256'],tested_steps=tested,
        first_step=first,physical_start_s=4e-7,physical_end_s=b.time,maximum_restart_difference=differences,
        wall_seconds=time.perf_counter()-start,
        note='Equivalence to stored accepted trajectory, not independent physics reference; no files modified in input case.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh output required')
    try:r=check(a.case);code=0
    except Exception as e:r=dict(status='FAIL_CHECKPOINT_READBACK',error=str(e),production_authorized=False);code=2
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(r,allow_nan=False));return code


if __name__=='__main__':raise SystemExit(main())
