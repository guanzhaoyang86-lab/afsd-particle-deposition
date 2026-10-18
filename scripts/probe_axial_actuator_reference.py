"""Two-kick manufactured actuator reference using existing PySPH integration.

Not a material model, not AFSD. Literal analytical reference values are only
checked after integrating; they never enter the dynamics callback.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.integrator import HostLeapfrog
from afsd_backend.pysph.contact_kick import kick
from afsd_backend.pysph.axial_actuator import control,validate_interval
from probe_rotating_deformable_block import jsonable

SOURCES=['scripts/probe_axial_actuator_reference.py','src/afsd_backend/pysph/axial_actuator.py',
    'src/afsd_backend/pysph/integrator.py','src/afsd_backend/pysph/contact_kick.py',
    'src/afsd_backend/pysph/contact_cone.py','docs/AXIAL_FEED_CONTRACT.md',
    'scripts/probe_rotating_deformable_block.py','third_party/PySPH/pysph/sph/integrator.py',
    'third_party/PySPH/pysph/sph/integrator_step.py']


def benchmark(normal_load):
    if normal_load not in (0.,2.):raise ValueError('registered manufactured normal load 0 or2N')
    pa=get_particle_array(name='solid',x=[0.],y=[0.],z=[2.],u=[1.],v=[-2.],w=[-.5],m=[2.])
    for name in ('ax','ay','az','au','av','aw','arho','ae','e'):
        if name not in pa.properties:pa.add_property(name)
    traces=[];records=[]
    def callback(t,dt):
        c=control(np.array([2.]),t,0.,2.,.5,.5,.25)
        fn=normal_load if t<.25 else 0.
        n=[[0.,0.,1.]] if fn else np.empty((0,3));fns=[fn] if fn else []
        walls=np.zeros((len(fns),3));v=np.array([pa.u[0],pa.v[0],pa.w[0]])
        r=kick(v,2.,np.array([4.,-8.,6.]),n,fns,walls,0.,dt,
            free_axes=c['free_axes'][0],prescribed_velocity=c['prescribed_velocity'][0])
        pa.au[:],pa.av[:],pa.aw[:]=((r['velocity']-v)/dt)
        traces.append(dict(time_mid_s=t,phase=c['phase'],result=r,normal_force_N=fn))
    integrator=HostLeapfrog(pa,callback)
    for step in range(2):
        validate_interval(step*.25,.25,.25);integrator.advance(step*.25,.25)
        records.append(dict(step=step+1,time_s=(step+1)*.25,
            position_m=[pa.x[0],pa.y[0],pa.z[0]],velocity_m_s=[pa.u[0],pa.v[0],pa.w[0]]))
    # Constants from the analytical statement; not obtained from kick/helper.
    expected=[dict(velocity_m_s=[1.5,-3.,-.5],position_m=[.3125,-.625,1.875],
        constraint_impulse=[0.,0.,-1.5 if normal_load==0 else -2.],
        constraint_work_J=.75 if normal_load==0 else 1.,other_work_J=5.5,
        delta_kinetic_J=6.25,normal_work_J=0. if normal_load==0 else -.25),
        dict(velocity_m_s=[2.,-4.,.25],position_m=[.75,-1.5,1.84375],constraint_impulse=[0.,0.,0.],
             constraint_work_J=0.,other_work_J=8.5625,delta_kinetic_J=8.5625,normal_work_J=0.)]
    for i,ref in enumerate(expected):
        r=traces[i]['result']
        for field in ('velocity_m_s','position_m'):
            np.testing.assert_allclose(records[i][field],ref[field],rtol=0,atol=1e-12)
        for field in ('constraint_impulse','constraint_work_J','other_work_J','delta_kinetic_J'):
            np.testing.assert_allclose(r[field],ref[field],rtol=0,atol=1e-12)
        np.testing.assert_allclose(sum(np.sum(s['normal_work_J']) for s in r['intervals']),ref['normal_work_J'],rtol=0,atol=1e-12)
        if np.sum(r['friction_heat_J'])!=0.:raise AssertionError('zero-mu friction not zero')
    return dict(status='PASS_LITERAL_AXIAL_ACTUATOR_REFERENCE_NOT_AFSD',normal_load_first_kick_N=normal_load,
        initial_kinetic_J=5.25,initial_velocity_m_s=[1.,-2.,-.5],mass_kg=2.,dt_s=.25,
        external_force_N=[4.,-8.,6.],traces=traces,records=records,independent_reference=expected,
        reference_scope='2kg manufactured point; NOT aluminum parameters, no plastic/thermal fields')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('fresh output required')
    start=time.perf_counter()
    result=dict(classification='MANUFACTURED_AXIAL_ACTUATOR_REFERENCE_NOT_AFSD',production_authorized=False,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES})
    try:result.update(status='PASS',cases=[benchmark(0.),benchmark(2.)]);code=0
    except Exception as error:result.update(status='FAIL',error=str(error));code=2
    result['wall_seconds']=time.perf_counter()-start
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(jsonable(result),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}));return code


if __name__=='__main__':raise SystemExit(main())
