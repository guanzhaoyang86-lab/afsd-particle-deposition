"""Small independent boundary-closure diagnostic, never a simulation."""
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

STRESS=np.array([[2.,.4,-.3],[.4,3.,.2],[-.3,.2,1.]])


def bubble(points):
    # Coordinates divided by ell=1m; amplitude=1m/s. NOT prescribed dynamics.
    x,y,z=np.asarray(points).T
    v=np.zeros((len(x),3));L=np.zeros((len(x),3,3))
    v[:,0]=x*(1-x*x)*(1-y*y)*(1-z*z)
    L[:,0,0]=(1-3*x*x)*(1-y*y)*(1-z*z)
    L[:,0,1]=-2*y*x*(1-x*x)*(1-z*z)
    L[:,0,2]=-2*z*x*(1-x*x)*(1-y*y)
    return v,L


def gauss_reference():
    # Fixed exact 3-point rule, independent of SPH and nodal quadrature.
    a=np.sqrt(3./5.);q=[-a,0.,a];w=[5./9.,8./9.,5./9.]
    idx=np.array(list(itertools.product(range(3),repeat=3)))
    points=np.array(q)[idx];volumes=np.prod(np.array(w)[idx],axis=1)
    _,L=bubble(points)
    return float(np.sum(volumes*np.einsum('ij,nij->n',STRESS,L)))


def evaluate(n):
    dx=2./(n-1);index=np.array(list(itertools.product(range(n),repeat=3)))
    points=np.linspace(-1.,1.,n)[index]
    weights=dx*np.where((index==0)|(index==n-1),.5,1.)
    V=np.prod(weights,axis=1);boundary_axes=np.sum((index==0)|(index==n-1),axis=1)
    pa=get_particle_array(x=points[:,0],y=points[:,1],z=points[:,2],m=V,
        rho=np.ones(len(points)),h=np.full(len(points),1.3*dx))
    ext=np.zeros_like(points);areas={}
    for axis in range(3):
        area=np.prod(weights[:,[j for j in range(3) if j!=axis]],axis=1)
        for side in (-1.,1.):
            on=points[:,axis]==side
            areas[f'{axis}:{side}']=float(area[on].sum())
            ext[on]+=area[on,None]*(STRESS@(side*np.eye(3)[axis]))
    fields={};states={};fi=None
    for name in ('affine','bubble'):
        if name=='affine':
            v=points.copy();exactL=np.broadcast_to(np.eye(3),(len(points),3,3));exact=48.
        else:v,exactL=bubble(points);exact=0.
        pa.u[:],pa.v[:],pa.w[:]=v.T
        i,j,g,L,_,B=corrected_operators(pa)
        fi=internal_force(pa,np.broadcast_to(STRESS,L.shape),i,j,g,correction=B)
        wsph=float(np.sum(V*np.einsum('ij,nij->n',STRESS,L)))
        wquad=float(np.sum(V*np.einsum('ij,nij->n',STRESS,exactL)))
        wsurf=float(np.sum(ext*v));wforce=-float(np.sum(fi*v))
        parts=dict(gradient_reconstruction_W=wsph-wquad,
            volume_quadrature_W=wquad-exact,surface_quadrature_W=exact-wsurf)
        fields[name]=dict(sph_volume_power_W=wsph,analytic_gradient_quadrature_W=wquad,
            analytic_continuum_power_W=exact,surface_power_W=wsurf,
            negative_internal_force_power_W=wforce,adjoint_identity_residual_W=wsph-wforce,
            closure_defect_W=wsph-wsurf,decomposition=parts,
            decomposition_identity_residual_W=sum(parts.values())-(wsph-wsurf))
        states[name+'_virtual_velocity']=v.copy();states[name+'_sph_gradient']=L
    residual=fi+ext;scale=float(np.linalg.norm(ext));groups={}
    for codim,label in enumerate(('interior','face','edge','corner')):
        mask=boundary_axes==codim
        groups[label]=dict(count=int(mask.sum()),residual_norm_N=float(np.linalg.norm(residual[mask])),
            normalized_by_total_traction=float(np.linalg.norm(residual[mask])/scale))
    inside=np.flatnonzero(boundary_axes==0)
    strongest=inside[np.argsort(np.linalg.norm(residual[inside],axis=1))[-5:]][::-1]
    record=dict(n=n,particles=len(points),dx_m=dx,physical_bounds_m=[-1.,1.],
        volume_m3=float(V.sum()),face_areas_m2=areas,stress_Pa=STRESS.tolist(),
        gauss_bubble_reference_W=gauss_reference(),fields=fields,node_groups=groups,
        nodal_relative_residual=float(np.linalg.norm(residual)/scale),
        interior_external_force_max_N=float(np.max(np.abs(ext[inside]))),
        strongest_interior=[dict(index=int(k),position_m=points[k].tolist(),
            residual_force_N=residual[k].tolist()) for k in strongest])
    states.update(position_m=points,volume_m3=V,boundary_codimension=boundary_axes,
                  internal_force_N=fi,external_traction_force_N=ext,residual_force_N=residual)
    return record,states


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter();rows=[]
    for n in (5,9):
        row,state=evaluate(n);np.savez_compressed(out/f'diagnostic_n{n}.npz',**state);rows.append(row)
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted((ROOT/'src/afsd_backend/pysph').glob('*.py'))}
    data=dict(status='DIAGNOSTIC_COMPLETED_NOT_PATCH_PASS',classification='VIRTUAL_WORK_FIXTURE',
        rows=rows,wall_seconds=time.perf_counter()-start,source_hashes=sources,
        force_or_trajectory_fix='NOT_IMPLEMENTED',afsd='NOT_RUN',production_authorized=False)
    (out/'report.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,indent=2))


if __name__=='__main__':main()
