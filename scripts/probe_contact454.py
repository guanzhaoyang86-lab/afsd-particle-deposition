"""Read-only frozen-input replay; every invocation requires a new output directory."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import contact_kick as uut
from afsd_backend.pysph.kick_diagnostics import clean

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'outputs/coupling_20260930/event_dx05_cfl01/REJECTED_contact_debug.json'
    c=json.loads(source.read_text())['failed_local_context'];face=c['faces'][0]
    frozen=json.loads((ROOT/'outputs/coupling_20260930/new_candidate_RED.json').read_text())
    cone=frozen['failed_cone'];n=np.array(face['normal'])
    # Independent Decimal plane basis for the captured normal (nz exactly zero).
    with localcontext() as ctx:
        ctx.prec=70
        nx,ny=[Decimal(str(a)) for a in n[:2]]
        norm=(nx*nx+ny*ny).sqrt(); nx/=norm;ny/=norm
        x,y,z=[Decimal(str(a)) for a in cone['base'][0]]
        q=(-ny*x+nx*y,z)
        tangent=[-ny*q[0],nx*q[0],q[1]]
        lam=Decimal(str(cone['coefficient']))/Decimal(str(c['mass_kg']))
        raw=[-a/lam for a in q]; mag=sum(a*a for a in raw).sqrt()
        limit=Decimal(str(cone['limits_N'][0]))
        if mag>limit:raw=[a*limit/mag for a in raw]
        reference=np.array([float(-ny*raw[0]),float(nx*raw[0]),float(raw[1])])
        tangent_reference=np.array([float(a) for a in tangent])
    projector=np.array(cone['projectors'])
    projected=projector[0]@np.array(cone['base'][0])
    normal_leakage=float(n@projected)
    report=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        uut_sha256=hashlib.sha256(Path(uut.__file__).read_bytes()).hexdigest(),
        cone_decimal_reference_N=reference,projected_velocity=projected,
        decimal_tangent_velocity=tangent_reference,projector_normal_leakage_m_s=normal_leakage)
    begin=time.perf_counter()
    fn=2.63e12*.0005**2*face['penetration']
    try:
        result=uut.kick(c['velocity_before'],c['mass_kg'],c['internal_force'],
            [n],[fn],[face['wall_velocity']],.3,c['dt_s'],free_axes=c['free_axes'],
            prescribed_velocity=c['prescribed_velocity'])
        report.update(status='COMPLETED_LOCAL_KICK',local_result=result)
    except uut.ContactConvergenceError as error:
        report.update(status='FAIL',error=str(error),failed_cone=getattr(error,'cone_context',None))
    report['wall_seconds']=time.perf_counter()-begin
    (out/'report.json').write_text(json.dumps(clean(report),indent=2)+'\n')
    print(json.dumps(clean({k:v for k,v in report.items() if k not in ('local_result','failed_cone')}),indent=2))
    return 1 if report['status']=='FAIL' else 0
if __name__=='__main__':raise SystemExit(main())
