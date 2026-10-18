"""Validate references/diagnostic identities, not a claim that patch passes."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from probe_boundary_closure import bubble,gauss_reference,evaluate


class BoundaryClosureDiagnosis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.cases=[evaluate(n) for n in (5,9)]

    def test_independent_bubble_reference_and_derivative(self):
        self.assertLess(abs(gauss_reference()),1e-12)
        x=np.array([[-.4,.2,.7],[.3,-.6,.1]])
        _,exact=bubble(x)
        for eps in (1e-4,1e-5,1e-6):
            for axis in range(3):
                offset=np.eye(3)[axis]*eps
                fd=(bubble(x+offset)[0]-bubble(x-offset)[0])/(2*eps)
                np.testing.assert_allclose(fd,exact[:,:,axis],atol=1e-8,rtol=0.)

    def test_fixed_geometry_traction_and_zero_boundary_virtual_velocity(self):
        for row,state in self.cases:
            with self.subTest(n=row['n']):
                self.assertAlmostEqual(row['volume_m3'],8.,delta=1e-12)
                for area in row['face_areas_m2'].values():self.assertAlmostEqual(area,4.,delta=1e-12)
                self.assertEqual(row['interior_external_force_max_N'],0.)
                np.testing.assert_array_equal(state['bubble_virtual_velocity'][state['boundary_codimension']>0],0.)
                np.testing.assert_allclose(state['external_traction_force_N'].sum(axis=0),0.,atol=1e-12)

    def test_affine_and_decomposition_identities(self):
        for row,_ in self.cases:
            self.assertAlmostEqual(row['fields']['affine']['sph_volume_power_W'],48.,delta=1e-10)
            self.assertAlmostEqual(row['fields']['affine']['surface_power_W'],48.,delta=1e-10)
            for values in row['fields'].values():
                self.assertLess(abs(values['adjoint_identity_residual_W']),1e-10)
                self.assertLess(abs(values['decomposition_identity_residual_W']),1e-10)

    def test_obstruction_is_observed_not_a_green_patch_test(self):
        for row,_ in self.cases:
            self.assertGreater(row['node_groups']['interior']['normalized_by_total_traction'],1e-8)
            self.assertGreater(abs(row['fields']['bubble']['closure_defect_W']),1e-8)
            self.assertEqual(row['fields']['bubble']['surface_power_W'],0.)


if __name__=='__main__':unittest.main(verbosity=2)
