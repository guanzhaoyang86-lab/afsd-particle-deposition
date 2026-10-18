"""Independent hand algebra / continuous ODE references for BE heat only."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.implicit_thermal import advance_implicit_thermal, pair_conductance


def arguments(n=2):
    x=np.column_stack((np.arange(n),np.zeros(n),np.zeros(n)))
    return dict(mass=np.ones(n),rho=np.ones(n),T=np.linspace(300.,400.,n),conductivity=.5,
        x=x,i=np.arange(n-1,dtype=int),j=np.arange(1,n,dtype=int),
        gradient=np.tile([1.,0.,0.],(n-1,1)),dt=1.,cp=1.,
        plastic_work_J=np.zeros(n),friction_raw_J=np.zeros(n),friction_rounding_bound_J=np.zeros(n),
        wall_row_W_K=np.zeros(n),wall_weighted_temperature_W=np.zeros(n),
        chi=.9,partition=.25,room=293.15,melt=855.15)


class ImplicitHeatTests(unittest.TestCase):
    def test_two_body_literal_inverse_exceeds_old_explicit_bound(self):
        r=advance_implicit_thermal(**arguments())
        np.testing.assert_allclose(r['temperature'],[1000/3,1100/3],atol=2e-12,rtol=0)
        self.assertEqual(r['former_explicit_dt_bound_s'],.25)
        self.assertAlmostEqual(r['ledger']['measured_thermal_increment_J'],0.,places=12)
        self.assertAlmostEqual(r['ledger']['internal_conduction_J'],0.,places=12)

    def test_same_end_time_against_independent_exponential(self):
        errors=[];exact=np.array([350-50*np.exp(-2),350+50*np.exp(-2)])
        for n in (10,20,40):
            a=arguments();a['dt']=1/n
            for step in range(n):a['T']=advance_implicit_thermal(**a)['temperature']
            errors.append(float(np.max(np.abs(a['T']-exact))))
        print('independent_exponential_same_end_errors_K',errors,flush=True)
        for ratio in np.array(errors[:-1])/errors[1:]:self.assertTrue(1.7<ratio<2.2)

    def test_one_body_wall_plus_source_literal(self):
        a=arguments(1);a.update(cp=2.,dt=.5,plastic_work_J=np.array([4.]),chi=1.,
            wall_row_W_K=np.array([3.]),wall_weighted_temperature_W=np.array([1500.]))
        r=advance_implicit_thermal(**a)
        np.testing.assert_allclose(r['temperature'],[1354/3.5],atol=2e-12,rtol=0)
        self.assertAlmostEqual(r['ledger']['contact_in_J'],.5*3*(500-1354/3.5),places=11)

    def test_source_partition_not_double_counted(self):
        a=arguments(1);a.update(cp=2.,plastic_work_J=np.array([10.]),friction_raw_J=np.array([4.]))
        r=advance_implicit_thermal(**a)
        np.testing.assert_array_equal(r['temperature'],[305.])
        self.assertAlmostEqual(r['ledger']['plastic_unthermalized_J'],1.)
        self.assertEqual(r['ledger']['friction_wall_heat_J'],3.)

    def test_signed_roundoff_friction_not_clipped(self):
        a=arguments(1);a['friction_raw_J'][0]=-1e-26;a['friction_rounding_bound_J'][0]=2e-26
        r=advance_implicit_thermal(**a)
        self.assertEqual(r['ledger']['friction_raw_J'],-1e-26)
        self.assertEqual(r['ledger']['friction_material_heat_J'],-.25e-26)
        a['friction_rounding_bound_J'][0]=0
        with self.assertRaisesRegex(ValueError,'not clipped'):advance_implicit_thermal(**a)

    def test_unequal_capacity_literal_conservation(self):
        a=arguments();a['cp']=np.array([2.,1.]);r=advance_implicit_thermal(**a)
        # [3,-1;-1,2] [T1,T2] = [600,400], determinant=5.
        np.testing.assert_allclose(r['temperature'],[320.,360.],atol=2e-12,rtol=0)
        self.assertAlmostEqual(r['ledger']['measured_thermal_increment_J'],0.,places=12)

    def test_pair_coefficient_independent_mixed_material_algebra(self):
        g=pair_conductance([2.,3.],[4.,5.],[6.,9.],[[0.,0.,0.],[2.,0.,0.]],
            np.array([0]),np.array([1]),[[.5,0.,0.]])
        # (6/20)*(216/15)*(1/4)=1.08 W/K.
        np.testing.assert_allclose(g,[1.08],rtol=1e-13,atol=0)

    def test_maximum_principle_for_large_dt_without_clipping(self):
        a=arguments(3);a['T']=np.array([400.,300.,500.]);a['dt']=100.
        r=advance_implicit_thermal(**a)
        self.assertTrue(np.all((r['temperature']>=300)&(r['temperature']<=500)))
        self.assertAlmostEqual(r['temperature'].sum(),1200.,places=10)

    def test_permutation_and_purity(self):
        a=arguments(3);old=copy.deepcopy(a);r=advance_implicit_thermal(**a);perm=np.array([2,0,1]);inv=np.argsort(perm)
        b=copy.deepcopy(a)
        for key in ('mass','rho','T','x','plastic_work_J','friction_raw_J','friction_rounding_bound_J',
                    'wall_row_W_K','wall_weighted_temperature_W'):b[key]=b[key][perm]
        b['i']=inv[a['i']];b['j']=inv[a['j']]
        np.testing.assert_allclose(advance_implicit_thermal(**b)['temperature'][inv],r['temperature'],atol=2e-12,rtol=0)
        for key in a:np.testing.assert_array_equal(a[key],old[key])

    def test_duplicate_or_bad_sign_or_nonradial_edges_rejected(self):
        for kind in ('duplicate','sign','nonradial','self','outside'):
            a=arguments()
            if kind=='duplicate':a.update(i=np.array([0,1]),j=np.array([1,0]),gradient=np.array([[1.,0,0],[-1.,0,0]]))
            elif kind=='sign':a['gradient']*=-1
            elif kind=='nonradial':a['gradient'][0,1]=.1
            elif kind=='self':a['j'][0]=0
            else:a['j'][0]=2
            with self.subTest(kind=kind),self.assertRaises(ValueError):advance_implicit_thermal(**a)

    def test_domain_invalid_parameters_and_no_silent_heat(self):
        changes=[('cp',0.),('dt',0.),('T',np.array([900.,400.])),('chi',2.),
            ('plastic_work_J',np.array([-1.,0.])),('wall_row_W_K',np.array([-1.,0.])),
            ('wall_weighted_temperature_W',np.array([1.,0.])),('mass',np.array([np.nan,1.]))]
        for key,value in changes:
            a=arguments();a[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):advance_implicit_thermal(**a)
        a=arguments(1);a['plastic_work_J'][0]=1000
        with self.assertRaisesRegex(ValueError,'NOT clipped'):advance_implicit_thermal(**a)


if __name__=='__main__':unittest.main(verbosity=2)
