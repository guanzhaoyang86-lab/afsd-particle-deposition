"""Bounded two-step restart check at the actual hold/lift transition.

This is restart equivalence, not an independent continuum solution. Input run
files are read-only. The original material/NNPS/contact path actually executes.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from probe_normal_deformable_block import Block,jsonable
from inspect_normal_contact import provenance
from inspect_rotating_deformable_block import load,file_hash,require,close


def check(directory):
    start=time.perf_counter();directory=Path(directory)
    m=provenance(directory,'COMPLETE_NORMAL_DEFORMABLE_COMPONENT_NOT_AFSD')
    history=load(directory/'history.json');checkpoint=directory/'prelift_accepted.npz'
    block=Block.from_checkpoint(checkpoint)
    require(block.step in (40,80) and block.time==2e-7,'hold/lift checkpoint time')
    start_step=block.step;maxima={};checks=[]
    for _ in range(3):
        r=jsonable(block.record());ref=history[block.step]
        require(r['particle_id']==ref['particle_id'],'persistent checkpoint IDs')
        require(r['tracer_type']==ref['tracer_type'] and r['as_rs_tag']==ref['as_rs_tag'],'persistent labels')
        for field,tol in [('position_m',1e-12),('velocity_m_s',1e-12),('density_kg_m3',1e-12),
            ('deviatoric_stress_Pa',1e-6),('ep',1e-12)]:
            defect=close(r[field],ref[field],tol,'restarted '+field)
            maxima[field]=max(maxima.get(field,0.),defect)
        for field in ref['cumulative']:
            close(r['cumulative'][field],ref['cumulative'][field],1e-15,'restarted ledger '+field)
        close(r['normal_contact']['gap_m'],ref['normal_contact']['gap_m'],1e-17,'restart normal gap')
        close(r['normal_contact']['normal_force_N'],ref['normal_contact']['normal_force_N'],1e-12,'restart normal force')
        close(r['tool_translation_velocity_m_s'],[0.,0.,20.],0.,'right-hand lift program')
        checks.append(block.step)
        if block.step<start_step+2: block.advance()
    return dict(status='PASS_NONZERO_PRELIFT_RESTART_TWO_STEPS_NOT_AFSD',production_authorized=False,
        source_sha256=m['source_sha256'],checkpoint_sha256=file_hash(checkpoint),
        manifest_sha256=file_hash(directory/'manifest.json'),script_sha256=file_hash(Path(__file__)),
        tested_steps=checks,physical_start_s=2e-7,physical_end_s=block.time,
        maximum_restart_difference=maxima,wall_seconds=time.perf_counter()-start,
        note='Stored trajectory as restart reference only; not a new physical benchmark or extended AFSD run.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh output required')
    try:r=check(a.case);code=0
    except Exception as e:r=dict(status='FAIL_RESTART',error=str(e),production_authorized=False);code=2
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(r,allow_nan=False));return code


if __name__=='__main__':raise SystemExit(main())
