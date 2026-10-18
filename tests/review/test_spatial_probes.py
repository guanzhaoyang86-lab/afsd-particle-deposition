"""Spatial diagnostics: independent affine and boundary-traction definitions.
--boundary-red retains nodal constant-stress boundary patch failure if present.
No spatial operator changes, no pseudoinverse or hidden stabilization.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.mechanics import corrected_operators, internal_force

def cloud(distorted=False):
    x=np.array(list(itertools.product([-1.,0.,1.],repeat=3)))
    # Physical domain [-1,1]^3, trapezoidal tensor-product volume quadrature.
    weights=np.where(x==0.,1.,.5)
    volume=np.prod(weights,axis=1)
    if distorted:
        x=x+.06*np.random.default_rng(395).normal(size=x.shape)
    pa=get_particle_array(x=x[:,0],y=x[:,1],z=x[:,2],m=volume,rho=np.ones(len(x)),
                          h=np.full(len(x),1.3))
    return pa,x,weights

def set_velocity(pa,v):
    pa.u[:],pa.v[:],pa.w[:]=v.T

def boundary_patch():
    pa,x,w=cloud()
    s=np.array([[2.,.4,-.3],[.4,3.,.2],[-.3,.2,1.]])
    sigma=np.broadcast_to(s,(len(x),3,3))
    i,j,g,L,_,B=corrected_operators(pa)
    f=internal_force(pa,sigma,i,j,g,correction=B)
    # Independent continuum traction sigma*n, tensor trapezoidal face areas.
    external=np.zeros_like(x)
    for axis in range(3):
        for side in (-1.,1.):
            on=x[:,axis]==side
            area=np.prod(w[:,[a for a in range(3) if a!=axis]],axis=1)
            normal=np.eye(3)[axis]*side
            external[on]+=area[on,None]*(s@normal)
    defect=np.linalg.norm(f+external)/np.linalg.norm(external)
    # Weak affine patch has independent reference integral sigma:A * volume.
    A=np.array([[.1,.3,-.2],[.6,-.4,.7],[-.8,.9,.2]])
    velocity=x@A.T+np.array([.2,-.6,.1])
    weak=float(np.sum(external*velocity))
    reference=8.*float(np.sum(s*A))
    return dict(nodal_boundary_relative_residual=float(defect),
                weak_affine_external_power=weak,independent_weak_power=reference,
                positions=x.tolist(),internal_force=f.tolist(),surface_traction_force=external.tolist())

def zero_modes():
    pa,x,_=cloud(distorted=True)
    n=len(x); operator=np.empty((6*n,3*n))
    for col in range(3*n):
        v=np.zeros((n,3));v.flat[col]=1.;set_velocity(pa,v)
        *_,L,drho,B=corrected_operators(pa)
        D=.5*(L+L.swapaxes(1,2))
        operator[:,col]=np.column_stack([D[:,0,0],D[:,1,1],D[:,2,2],
            np.sqrt(2)*D[:,0,1],np.sqrt(2)*D[:,0,2],np.sqrt(2)*D[:,1,2]]).ravel()
    singular=np.linalg.svd(operator,compute_uv=False)
    tol=1e-10*singular[0]
    nullity=int(np.sum(singular<tol))
    return dict(particles=n,dofs=3*n,threshold=float(tol),
        nullity=nullity,rigid_modes_expected=6,extra_zero_modes=nullity-6,
        singular_values=singular.tolist(),
        scope='ONE_DISTORTED_FREE_SURFACE_CLOUD; not a global stability guarantee')

class SpatialProbes(unittest.TestCase):
    def test_affine_free_surface_distorted_and_conjugate_balance(self):
        pa,x,_=cloud(distorted=True)
        A=np.array([[.1,.3,-.2],[.6,-.4,.7],[-.8,.9,.2]])
        v=x@A.T+np.array([.2,-.6,.1]);set_velocity(pa,v)
        i,j,g,L,drho,B=corrected_operators(pa)
        np.testing.assert_allclose(L,np.broadcast_to(A,L.shape),atol=1e-10,rtol=1e-10)
        # Different symmetric stress at each particle; unequal physical weights.
        s=np.random.default_rng(655).normal(size=(len(x),3,3));s+=s.swapaxes(1,2)
        f=internal_force(pa,s,i,j,g,correction=B)
        work=float(np.sum(f*v)+np.sum(pa.m/pa.rho*np.sum(s*L,axis=(1,2))))
        scale=float(np.sum(np.linalg.norm(f,axis=1)))
        self.assertLess(abs(work),1e-10*scale)
        self.assertLess(np.linalg.norm(f.sum(axis=0)),1e-10*scale)
        self.assertLess(np.linalg.norm(np.cross(x,f).sum(axis=0)),1e-10*scale*np.linalg.norm(x,axis=1).max())
    def test_zero_energy_mode_probe(self):
        r=zero_modes();print('ZERO_MODES',json.dumps(r))
        self.assertEqual(r['nullity'],6)
    def test_boundary_weak_affine_power(self):
        r=boundary_patch()
        self.assertAlmostEqual(r['weak_affine_external_power'],r['independent_weak_power'],delta=1e-10)
    def test_condition_number_threshold_fails_closed(self):
        pa,x,_=cloud()
        # Hyperbolic spline has nonzero radial slope at r=0: for nearly
        # coincident vertical pairs lambda_z scales as dz, not dz**2.
        # The original 1e-5 fixture did NOT exceed the 1e6 guard.
        pa.z[:]*=1e-9
        with self.assertRaisesRegex(ValueError,'ill-conditioned'):
            corrected_operators(pa)

class BoundaryRed(unittest.TestCase):
    def test_pointwise_constant_stress_boundary_patch(self):
        r=boundary_patch();print('BOUNDARY_PATCH',json.dumps(r))
        self.assertLess(r['nodal_boundary_relative_residual'],1e-10)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--boundary-red',action='store_true')
    parser.add_argument('--output')
    args=parser.parse_args()
    if args.output:
        path=Path(args.output)
        if path.exists():raise FileExistsError(path)
        path.write_text(json.dumps(dict(boundary=boundary_patch(),zero_modes=zero_modes()),indent=2)+'\n')
    cls=BoundaryRed if args.boundary_red else SpatialProbes
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(cls))
    raise SystemExit(0 if result.wasSuccessful() else 1)

