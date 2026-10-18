"""Literal independent-reference invariants; no tested HHT/backend imports."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from reference_hht_block import fixture,wall_program,pair_fields,boundary_fields,packed_initial,rhs,unpack


class HHTReference(unittest.TestCase):
    def test_literal_mass_geometry_and_initial_state(self):
        d=fixture()
        self.assertEqual(d['position_m'].shape,(27,3))
        self.assertAlmostEqual(float(d['mass_kg'].sum()),.0027,delta=1e-16)
        self.assertEqual([len(ids) for ids in d['face_nodes']],[9,9])
        np.testing.assert_array_equal(d['velocity_m_s'][:,2],.02)
        self.assertTrue(np.all(d['deviatoric_stress_Pa'][:,0,0]==20000))
        self.assertEqual(float(d['temperature_K'].min()),400.)
        self.assertEqual(float(d['temperature_K'].max()),402.)

    def test_literal_initial_contact_and_reversal(self):
        d=fixture(); p=d['parameters']; x=d['position_m']; v=d['velocity_m_s']; T=d['temperature_K']
        a=boundary_fields(0.,x,v,T,d,False)
        np.testing.assert_allclose(a['face_force'],[[126.,0.,25.2],[-126.,0.,25.2]],rtol=0,atol=5e-11)
        np.testing.assert_allclose(a['friction'].sum(axis=1),[4.536,4.536],rtol=0,atol=1e-11)
        np.testing.assert_allclose(a['wall_heat'].sum(axis=1),[2.25,4.05],rtol=0,atol=1e-13)
        left,wl=wall_program(1e-6,p,False); right,wr=wall_program(1e-6,p,True)
        np.testing.assert_array_equal(left,right)
        self.assertEqual(wl[0,2],.2); self.assertEqual(wr[0,2],-.2)
        before=a['normal_material_work']+a['tangent_material_work']
        np.testing.assert_allclose(before-a['wall_work']+a['friction'].sum(axis=1),0.,atol=1e-12)

    def test_independent_current_pair_power_and_affine(self):
        d=fixture(); x=d['position_m']; A=np.array([[.3,.1,-.2],[-.2,.4,.1],[.1,.2,-.3]])
        v=x@A.T+np.array([.1,-.2,.3]); S=d['deviatoric_stress_Pa']; rho=d['density_kg_m3']
        f=pair_fields(x,v,S,rho,d)
        np.testing.assert_allclose(f['gradient'],np.broadcast_to(A,(27,3,3)),rtol=0,atol=1e-13)
        self.assertEqual(len(f['a']),305)
        np.testing.assert_allclose(f['force'].sum(axis=0),0.,atol=1e-13)
        np.testing.assert_allclose(np.cross(x,f['force']).sum(axis=0),0.,atol=1e-14)
        defect=np.sum(f['force']*v)+np.sum((d['mass_kg']/rho)[:,None,None]*S*f['gradient'])
        self.assertLess(abs(defect),1e-12)
        self.assertTrue(np.all(f['conductance']>=0))

    def test_rhs_preserves_nonzero_initial_history_and_heat(self):
        d=fixture(); y=packed_initial(d); before=y.copy(); derivative=rhs(0.,y,d,False)
        np.testing.assert_array_equal(y,before)
        x,v,S,rho,T,ell=unpack(derivative)
        np.testing.assert_allclose(x,d['velocity_m_s'],atol=0,rtol=0)
        self.assertGreater(float(np.max(np.abs(v))),1e4)
        self.assertGreater(float(np.max(T)),0.)
        self.assertAlmostEqual(float(np.sum(d['mass_kg']*900*T)),6.3+.5*9.072,delta=1e-10)


if __name__=='__main__': unittest.main(verbosity=2)
