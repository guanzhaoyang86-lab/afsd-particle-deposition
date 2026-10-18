"""Independent SI references for a smooth frozen-rule contact Jacobian.

Expected forces, dense kinematics and complete DeltaJ residuals below do not
call production helpers. No material update or physical trajectory is run.
Pre-registration: docs/PROCESS_RESPONSE_TEST_CONTRACT.md.
"""
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from afsd_backend.pysph import contact_path_response as tested

ATOL, RTOL = 2e-8, 2e-6
EXACT_ATOL, EXACT_RTOL = 2e-12, 2e-12
EPSILONS = (2e-5, 1e-5)
ROTATION = np.array([[0., -.8, .6], [1., 0., 0.], [0., .6, .8]])


def emit(name, **values):
    print('CONTACT_RESPONSE_EVIDENCE ' + json.dumps(dict(name=name,
          classification='NUMERICAL_COMPONENT_NOT_TRAJECTORY', **values),
          sort_keys=True, allow_nan=False), flush=True)


def serialized(value):
    return json.dumps(value, sort_keys=True,
                      default=lambda a: np.asarray(a).tolist(), allow_nan=False)


def literal_force(normal, slip, fn, mu, sign):
    """Only algebraic force values, not a derivative or tested routine."""
    return fn*normal - sign*mu*fn*slip/np.linalg.norm(slip)


def central_columns(function, value, epsilon):
    eye = np.eye(3)
    return np.column_stack([(function(value+epsilon*e)-function(value-epsilon*e))
                            / (2*epsilon) for e in eye])


def manufactured_path(dt, force_weight, sign):
    return dict(dt=dt, force_weight=force_weight, sign=sign,
                n=ROTATION @ np.array([0., 0., 1.]),
                x0=ROTATION @ np.array([.03, -.07, -.02]),
                v0=ROTATION @ np.array([.9, .5, .03]),
                acc0=ROTATION @ np.array([.4, -.2, .06]),
                acc1=ROTATION @ np.array([-.6, .5, -.03]),
                wall0=ROTATION @ np.array([.2, -.3, .01]),
                wall_change=ROTATION @ np.array([.1, .07, 0.]),
                fn0=6.2, kn=31., mu=.37)


def literal_row(case, acceleration, fraction):
    h=case['dt']; a=case['force_weight']; gamma=1.5-a
    beta=(2-a)**2/4; d=beta-gamma/2; s=fraction
    da=acceleration-case['acc0']
    dv=h*((1-gamma)*case['acc0']+gamma*acceleration)
    velocity=case['v0']+s*dv+6*h*d*da*s*(1-s)
    position=case['x0']+h*(s*case['v0']+.5*s*s*dv+h*d*da*(3*s*s-2*s**3))
    normal=case['n']; relative=velocity-case['wall0']-s*case['wall_change']
    slip=relative-normal*(normal @ relative)
    fn=case['fn0']-case['kn']*(normal @ (position-case['x0']))
    if fn <= 0 or np.linalg.norm(slip) <= .1:
        raise AssertionError('reference path left its preregistered smooth domain')
    return dict(normal=normal.copy(), tangential_relative_velocity_m_s=slip,
                normal_force_N=fn, normal_stiffness_N_m=case['kn'],
                friction=case['mu'], branch_sign=case['sign'])


def force_from_row(row):
    return literal_force(row['normal'], row['tangential_relative_velocity_m_s'],
                         row['normal_force_N'], row['friction'], row['branch_sign'])


def quadrature_rows(case, acceleration, order):
    abscissa, weights=np.polynomial.legendre.leggauss(order)
    return [dict(literal_row(case, acceleration, .5*(s+1)),
                 fraction=.5*(s+1), weight=.5*w) for s,w in zip(abscissa,weights)]


def literal_integrals_and_root(case, acceleration, order=64):
    h=case['dt']; a=case['force_weight']; gamma=1.5-a
    rows=quadrature_rows(case, acceleration, order)
    impulse=h*sum(row['weight']*force_from_row(row) for row in rows)
    f0=force_from_row(literal_row(case, acceleration, 0.))
    f1=force_from_row(literal_row(case, acceleration, 1.))
    delta_impulse=impulse-.5*h*(f0+f1)
    # Crucially this reference differentiates the complete, independently
    # evaluated residual; it does not reuse the tested Jacobian formula.
    root=-a*f1-(1-a)*f0-delta_impulse/(gamma*h)
    return impulse, f1, root


class ContactPathResponseReferences(unittest.TestCase):
    def exact(self, value, expected):
        np.testing.assert_allclose(value, expected, atol=EXACT_ATOL, rtol=EXACT_RTOL)

    def finite_difference(self, value, expected):
        np.testing.assert_allclose(value, expected, atol=ATOL, rtol=RTOL)

    def test_01_literal_point_matrices_both_signs_and_zero_friction(self):
        n=np.array([0.,0.,1.]);u=np.array([3.,4.,0.])
        for sign in (1,-1):
            with self.subTest(sign=sign):
                result=tested.point_response(n,u,5.,7.,.2,sign)
                expected_x=np.array([[0.,0.,.84*sign],[0.,0.,1.12*sign],[0.,0.,-7.]])
                expected_v=sign*np.array([[-.128,.096,0.],[.096,-.072,0.],[0.,0.,0.]])
                self.exact(result['Kx_N_m'],expected_x)
                self.exact(result['Kv_N_s_m'],expected_v)
                self.assertFalse(result['accepted']);self.assertFalse(result['restart_allowed'])
                self.assertFalse(result['continuous_path_differentiability_certified'])
                self.exact(result['Kv_N_s_m'] @ u,np.zeros(3))
                self.exact(result['Kv_N_s_m'] @ n,np.zeros(3))
        zero=tested.point_response(n,u,5.,7.,0.,1)
        self.exact(zero['Kx_N_m'],np.diag([0.,0.,-7.]))
        self.exact(zero['Kv_N_s_m'],np.zeros((3,3)))

    def test_02_point_force_value_finite_differences_include_normal_gap(self):
        worst_x=0.;worst_v=0.
        for rotation in (np.eye(3),ROTATION):
            n=rotation @ np.array([0.,0.,1.]);u=rotation @ np.array([3.,4.,0.])
            for sign in (1,-1):
                result=tested.point_response(n,u,5.,7.,.2,sign)
                def force_x(x):return literal_force(n,u,5.-7.*(n@x),.2,sign)
                def force_v(v):return literal_force(n,v-n*(n@v),5.,.2,sign)
                for epsilon in EPSILONS:
                    dx=central_columns(force_x,np.zeros(3),epsilon)
                    dv=central_columns(force_v,u,epsilon)
                    self.finite_difference(result['Kx_N_m'],dx)
                    self.finite_difference(result['Kv_N_s_m'],dv)
                    worst_x=max(worst_x,float(np.max(abs(result['Kx_N_m']-dx))))
                    worst_v=max(worst_v,float(np.max(abs(result['Kv_N_s_m']-dv))))
        emit('point_actual_force_central_differences',maximum_Kx_error_N_m=worst_x,
             maximum_Kv_error_N_s_m=worst_v,epsilons=list(EPSILONS))

    def test_03_rotation_covariance_and_signed_transverse_sensitivity(self):
        n=np.array([0.,0.,1.]);u=np.array([3.,4.,0.])
        self.exact(ROTATION.T @ ROTATION,np.eye(3));self.exact(np.linalg.det(ROTATION),1.)
        results={}
        for sign in (1,-1):
            plain=tested.point_response(n,u,5.,7.,.2,sign)
            rotated=tested.point_response(ROTATION@n,ROTATION@u,5.,7.,.2,sign)
            for key in ('Kx_N_m','Kv_N_s_m'):
                self.exact(rotated[key],ROTATION @ plain[key] @ ROTATION.T)
            results[sign]=plain
        tangent=np.array([-.8,.6,0.])
        self.assertLess(tangent @ results[1]['Kv_N_s_m'] @ tangent,0.)
        self.assertGreater(tangent @ results[-1]['Kv_N_s_m'] @ tangent,0.)
        separation=float(np.max(abs(results[1]['Kv_N_s_m']-results[-1]['Kv_N_s_m'])))
        self.assertGreater(separation,1e-3)
        emit('private_wrong_sign_sensitivity',maximum_Kv_separation_N_s_m=separation)

    def test_04_complete_path_impulse_endpoint_and_root_finite_differences(self):
        maxima=np.zeros(3);quad_error=0.;coarse_fine=0.;endpoint_omission=0.
        samples=0;minimum_fn=float('inf');minimum_slip=float('inf')
        keys=('actual_impulse_acceleration_response_kg_s',
              'endpoint_force_acceleration_response_kg','root_contact_response_kg')
        for h in (.2,.7):
            for force_weight in (.9,1.,2./3.):
                for sign in (1,-1):
                    with self.subTest(h=h,force_weight=force_weight,sign=sign):
                        case=manufactured_path(h,force_weight,sign);acc=case['acc1']
                        nodes=quadrature_rows(case,acc,16);end=literal_row(case,acc,1.)
                        original=serialized((nodes,end))
                        result=tested.path_response(nodes,end,h,weight=force_weight)
                        self.assertEqual(serialized((nodes,end)),original)
                        self.assertFalse(result['accepted']);self.assertFalse(result['restart_allowed'])
                        self.assertFalse(result['physical_force_changed'])
                        self.assertFalse(result['continuous_path_differentiability_certified'])
                        self.assertFalse(result['quadrature_selection_derivative_included'])
                        gamma=1.5-force_weight
                        self.exact(result[keys[2]],(1/(2*gamma)-force_weight)*result[keys[1]]
                                   -result[keys[0]]/(gamma*h))
                        # Dropping the endpoint term of DeltaJ gives this wrong root.
                        missing=-force_weight*result[keys[1]]-result[keys[0]]/(gamma*h)
                        endpoint_omission=max(endpoint_omission,float(np.max(abs(missing-result[keys[2]]))))
                        center64=literal_integrals_and_root(case,acc,64)[0]
                        center128=literal_integrals_and_root(case,acc,128)[0]
                        np.testing.assert_allclose(center64,center128,atol=1e-12,rtol=1e-10)
                        quad_error=max(quad_error,float(np.max(abs(center64-center128))))
                        derivatives=[]
                        for epsilon in EPSILONS:
                            block=[]
                            for index,key in enumerate(keys):
                                fd=central_columns(lambda value:literal_integrals_and_root(case,value,64)[index],acc,epsilon)
                                self.finite_difference(result[key],fd);block.append(fd)
                                maxima[index]=max(maxima[index],float(np.max(abs(result[key]-fd))))
                            derivatives.append(block)
                        for low,high in zip(*derivatives):
                            self.finite_difference(low,high)
                            coarse_fine=max(coarse_fine,float(np.max(abs(low-high))))
                        for row in quadrature_rows(case,acc,128):
                            minimum_fn=min(minimum_fn,row['normal_force_N'])
                            minimum_slip=min(minimum_slip,float(np.linalg.norm(row['tangential_relative_velocity_m_s'])))
                        samples+=1
        self.assertEqual(samples,12);self.assertGreater(endpoint_omission,1e-3)
        emit('complete_independent_Gauss64_root_difference',cases=samples,
             maximum_impulse_response_error_kg_s=float(maxima[0]),
             maximum_endpoint_response_error_kg=float(maxima[1]),
             maximum_root_response_error_kg=float(maxima[2]),
             maximum_Gauss64_128_impulse_difference_Ns=quad_error,
             maximum_two_epsilon_matrix_difference=coarse_fine,
             missing_endpoint_correction_separation_kg=endpoint_omission,
             minimum_normal_force_N=minimum_fn,minimum_full_slip_m_s=minimum_slip)

    def test_05_path_rotation_covariance_and_detached_results(self):
        case=manufactured_path(.7,.9,-1);nodes=quadrature_rows(case,case['acc1'],16)
        end=literal_row(case,case['acc1'],1.);before=serialized((nodes,end))
        result=tested.path_response(nodes,end,.7)
        rotated_nodes=copy.deepcopy(nodes);rotated_end=copy.deepcopy(end)
        for row in rotated_nodes+[rotated_end]:
            for key in ('normal','tangential_relative_velocity_m_s'):row[key]=ROTATION @ row[key]
        rotated=tested.path_response(rotated_nodes,rotated_end,.7)
        keys=('actual_impulse_acceleration_response_kg_s',
              'endpoint_force_acceleration_response_kg','root_contact_response_kg')
        for key in keys:self.exact(rotated[key],ROTATION @ result[key] @ ROTATION.T)
        expected=result['root_contact_response_kg'].copy()
        for key in keys:result[key][...]=123.
        self.assertEqual(serialized((nodes,end)),before)
        repeated=tested.path_response(nodes,end,.7)
        self.exact(repeated['root_contact_response_kg'],expected)

    def test_06_point_invalid_input_rejection_without_regularization(self):
        valid=dict(normal=np.array([0.,0.,1.]),slip_m_s=np.array([3.,4.,0.]),
                   normal_force_N=5.,normal_stiffness_N_m=7.,friction=.2,branch_sign=1)
        changes=[('normal',[0.,0.,2.]),('normal',[0.,0.,np.nan]),('normal',[0.,1.]),
                 ('slip_m_s',[0.,0.,0.]),('slip_m_s',[1.,0.,.1]),('slip_m_s',[np.inf,0.,0.]),
                 ('normal_force_N',0.),('normal_force_N',-1.),('normal_force_N',np.inf),
                 ('normal_stiffness_N_m',0.),('normal_stiffness_N_m',-1.),
                 ('friction',-1.),('friction',np.nan),('branch_sign',0),('branch_sign',2),
                 ('branch_sign',True),('branch_sign',np.bool_(True))]
        before=serialized(valid)
        for key,value in changes:
            with self.subTest(key=key,value=str(value)):
                with self.assertRaises(ValueError):tested.point_response(**dict(valid,**{key:value}))
        self.assertEqual(serialized(valid),before)
        tiny=tested.point_response([0.,0.,1.],[1e-10,0.,0.],5.,7.,.2)
        self.exact(tiny['Kv_N_s_m'][1,1],-1e10)

    def test_07_path_invalid_quadrature_and_fixed_geometry_rejection(self):
        case=manufactured_path(.2,.9,1);nodes=quadrature_rows(case,case['acc1'],16)
        end=literal_row(case,case['acc1'],1.);before=serialized((nodes,end))
        for dt in (0.,-1.,np.inf):
            with self.subTest(dt=dt):
                with self.assertRaises(ValueError):tested.path_response(nodes,end,dt)
        for weight in (.5,1.1,np.nan):
            with self.subTest(weight=weight):
                with self.assertRaises(ValueError):tested.path_response(nodes,end,.2,weight=weight)
        with self.assertRaises(ValueError):tested.path_response([],end,.2)
        alterations=[('fraction',-.1),('fraction',1.1),('fraction',np.nan),
                     ('weight',0.),('weight',-1.),('weight',np.nan),('weight',.8),
                     ('normal',-case['n']),('normal_stiffness_N_m',32.),('friction',.38),
                     ('tangential_relative_velocity_m_s',np.zeros(3)),('branch_sign',0)]
        for key,value in alterations:
            altered=copy.deepcopy(nodes);altered[0][key]=value
            with self.subTest(key=key,value=str(value)):
                with self.assertRaises(ValueError):tested.path_response(altered,end,.2)
        self.assertEqual(serialized((nodes,end)),before)


if __name__=='__main__':
    sources=['src/afsd_backend/pysph/contact_path_response.py',
             'src/afsd_backend/pysph/hht_path.py',
             'tests/review/test_contact_path_response.py',
             'docs/PROCESS_RESPONSE_TEST_CONTRACT.md']
    emit('source_receipt_before_execution',sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                                 for name in sources})
    unittest.main(verbosity=2)
