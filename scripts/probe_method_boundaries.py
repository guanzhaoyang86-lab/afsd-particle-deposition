"""Fixed physical cube: independent traction/work probe, no operator edits."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.mechanics import corrected_operators,internal_force


def cube(n):
    spacing=2./(n-1)
    axes=np.linspace(-1.,1.,n)
    index=np.array(list(itertools.product(range(n),repeat=3)))
    x=axes[index]
    weights=np.where((index==0)|(index==n-1),.5,1.)*spacing
    volume=np.prod(weights,axis=1)
    pa=get_particle_array(x=x[:,0],y=x[:,1],z=x[:,2],m=volume,
        rho=np.ones(len(x)),h=np.full(len(x),1.3*spacing))
    S=np.array([[2.,.4,-.3],[.4,3.,.2],[-.3,.2,1.]])
    i,j,g,L,_,B=corrected_operators(pa)
    force=internal_force(pa,np.broadcast_to(S,(len(x),3,3)),i,j,g,correction=B)
    external=np.zeros_like(x)
    for axis in range(3):
        area=np.prod(weights[:,[a for a in range(3) if a!=axis]],axis=1)
        for sign in (-1.,1.):
            on=x[:,axis]==sign
            external[on]+=area[on,None]*(S@(sign*np.eye(3)[axis]))
    # Analytic continuum surface work = 8*trace(S)=48 for v=x AND v=x^3.
    # Cubic field has div-stress-work integral sum_i 3*Sii*int xi^2 dV.
    affine_internal=-float(np.sum(force*x))
    cubic_internal=-float(np.sum(force*x**3))
    return dict(n=n,particles=len(x),spacing_m=spacing,physical_cube_m=[-1.,1.],
        volume_m3=float(volume.sum()),stress_Pa=S.tolist(),
        nodal_relative_residual=float(np.linalg.norm(force+external)/np.linalg.norm(external)),
        affine_internal_work=affine_internal,affine_surface_work=float(np.sum(external*x)),
        cubic_internal_work=cubic_internal,cubic_surface_work=float(np.sum(external*x**3)),
        independent_affine_work=48.,independent_cubic_work=48.,
        cubic_relative_error=abs(cubic_internal-48.)/48.,
        net_internal_force=np.sum(force,axis=0).tolist(),
        net_internal_torque=np.sum(np.cross(x,force),axis=0).tolist())


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter();results=[cube(n) for n in (3,5,7)]
    sources={name:hashlib.sha256((ROOT/'src/afsd_backend/pysph'/name).read_bytes()).hexdigest()
             for name in ('material.py','mechanics.py','solver.py','geometry.py')}
    data=dict(status='DIAGNOSTIC_COMPLETED_NOT_METHOD_CERTIFIED',results=results,
        wall_seconds=time.perf_counter()-start,source_hashes=sources,
        warning='weak affine exactness does not imply nodal traction consistency; no thresholds relaxed',
        real_afsd='NOT_RUN',constitutive_change='NOT_IMPLEMENTED')
    (out/'report.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
    for a in results:
        if abs(a['affine_internal_work']-48.)>1e-10 or abs(a['volume_m3']-8.)>1e-12:
            raise AssertionError('fixed-domain/affine diagnostic failed')


if __name__=='__main__':main()
