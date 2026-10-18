"""Independent literal uneliminated-block oracles; never a physical map."""
from pathlib import Path
import copy,sys,unittest
from unittest import mock
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.material_stage_response import MaterialStageResponse

FIELDS=('position_m','velocity_m_s','density_kg_m3','temperature_K','mean_velocity_m_s')
def scalar():
    D=np.zeros((1,9,3));D[0,0,0]=1
    return dict(mass=np.array([4.]),volume=np.array([2.]),gradient_map=D,
        material_tangent=np.zeros((1,9,9)),material_temperature=np.zeros((1,9)),
        bulk=3.,rho_reference=5.,rho_parent=np.array([5.]),dt=2.,divergence=np.zeros(1),
        free_axes=np.ones((1,3),bool),particle_ids=np.array([17]))
def coupled():
    D=np.zeros((2,9,6));D[0,0,0]=1;D[0,4,3]=-1
    C=np.zeros((2,9,9));C[0,0,0]=2;C[0,0,4]=2;C[0,4,4]=4;C[0,8,0]=-2;C[0,8,4]=-6
    ST=np.zeros((2,9));ST[0,[0,4,8]]=[-4,-4,8]
    return dict(mass=np.array([2.,3.]),volume=np.ones(2),gradient_map=D,material_tangent=C,
        material_temperature=ST,bulk=2.,rho_reference=1.,rho_parent=np.array([.5,.5]),dt=2.,
        divergence=np.zeros(2),free_axes=np.ones((2,3),bool),particle_ids=np.array([17,29]))
def fields(n):
    return dict(position_m=np.zeros((n,3)),velocity_m_s=np.zeros((n,3)),
        density_kg_m3=np.full(n,5.),temperature_K=np.full(n,400.),mean_velocity_m_s=np.zeros((n,3)))
def coupled_fields():
    c=fields(2);g=copy.deepcopy(c);g['mean_velocity_m_s'][:,0]=[1,2]
    g['density_kg_m3'][0]+=.5;g['temperature_K'][0]+=.25
    return c,g

class MaterialStageResponseTests(unittest.TestCase):
    def close(self,a,b):np.testing.assert_allclose(a,b,rtol=2e-12,atol=2e-13)

    def test_01_pressure_sign_and_literal_scalar_schur(self):
        p=MaterialStageResponse(**scalar());c=fields(1);g=copy.deepcopy(c)
        g['mean_velocity_m_s'][0,0]=.5;g['density_kg_m3'][0]+=2
        got=p.apply(c,g);d=got['correction']
        self.close(d['mean_velocity_m_s'],[[.32,0,0]]);self.close(d['density_kg_m3'],[-1.2])
        self.close(d['position_m'],[[-.36,0,0]]);self.close(d['velocity_m_s'],[[-.36,0,0]])
        self.close(4*d['mean_velocity_m_s'][0,0]-.6*d['density_kg_m3'][0],2)
        self.assertLessEqual(got['diagnostics']['backward_error_inf'],1e-12)

    def test_02_exact_continuity_and_double_temperature_half(self):
        kw=scalar();kw.update(mass=np.array([2.]),volume=np.ones(1),bulk=2.,rho_reference=3.,
            rho_parent=np.array([3.]),dt=.5,divergence=np.array([.8]))
        kw['material_tangent'][0,0,0]=4;kw['material_temperature'][0,0]=-4
        p=MaterialStageResponse(**kw);c=fields(1);g=copy.deepcopy(c)
        g['mean_velocity_m_s'][0,0]=2/5;g['density_kg_m3'][0]+=3/5;g['temperature_K'][0]+=1/5
        d=p.apply(c,g)['correction']
        # Independent uneliminated block, not p.matrix or p.J*.
        ref=np.linalg.solve(np.array([[5/2,-1/12,-1/4],[25/24,1,0],[0,0,1.]]),[4/5,3/5,1/5])
        self.close(ref,[1296/3725,177/745,1/5])
        self.close([d['mean_velocity_m_s'][0,0],d['density_kg_m3'][0],d['temperature_K'][0]],ref)
        self.close(p.JT[:,0],[1,0,0]);self.close(p.Jr[:,0],[1/3,0,0]);self.close(p.kappa,[-25/24])

    def test_03_full_uneliminated_nonsymmetric_two_particle_reference(self):
        p=MaterialStageResponse(**coupled());c,g=coupled_fields();r=p.apply(c,g);d=r['correction']
        B=np.array([[3,-1,-1,-1],[0,5,1,1],[1,-1,1,0],[0,0,0,1.]])
        ref=np.linalg.solve(B,[2,6,.5,.25])
        self.close(ref,[27/22,95/88,31/88,.25])
        self.close([d['mean_velocity_m_s'][0,0],d['mean_velocity_m_s'][1,0],d['density_kg_m3'][0],d['temperature_K'][0]],ref)
        self.close(p.matrix[np.ix_([0,3],[0,3])],[[4,-2],[-1,6]])
        self.assertNotEqual(p.matrix[0,3],p.matrix[3,0])
        self.close(d['position_m'][:,0],2*(ref[:2]-[1,2]));self.close(d['velocity_m_s'][:,0],2*(ref[:2]-[1,2]))
        self.close(r['diagnostics']['algebraic_reaction_impulse_Ns'],np.zeros(6))

    def test_04_prescribed_column_density_feedback_and_reaction(self):
        kw=coupled();kw['free_axes'][1,0]=False;p=MaterialStageResponse(**kw);c,g=coupled_fields()
        r=p.apply(c,g);d=r['correction'];j=r['diagnostics']['algebraic_reaction_impulse_Ns']
        self.close(d['mean_velocity_m_s'][:,0],[27/16,2]);self.close(d['density_kg_m3'],[13/16,0])
        self.close(j[[0,3]],[0,81/16]);self.close(j[[1,2,4,5]],np.zeros(4))
        self.assertEqual(r['diagnostics']['reaction_scope'],'SEARCH_ALGEBRA_NOT_PHYSICAL_REACTION')

    def test_05_all_prescribed_empty_free_block(self):
        kw=coupled();kw['free_axes'][:]=False;c,g=coupled_fields()
        r=MaterialStageResponse(**kw).apply(c,g);d=r['correction']
        self.close(d['mean_velocity_m_s'][:,0],[1,2]);self.close(d['density_kg_m3'],[1.5,0])
        self.close(r['diagnostics']['algebraic_reaction_impulse_Ns'][[0,3]],[-11/4,23/4])
        self.assertEqual(r['diagnostics']['backward_error_inf'],0)
        self.close(d['position_m'],np.zeros((2,3)));self.close(d['velocity_m_s'],np.zeros((2,3)))

    def test_06_pure_detached_arrays_and_metadata(self):
        kw=coupled();before=copy.deepcopy(kw);p=MaterialStageResponse(**kw);c,g=coupled_fields();c0=copy.deepcopy(c);g0=copy.deepcopy(g)
        a=p.apply(c,g);saved=copy.deepcopy(a)
        for key in FIELDS:np.testing.assert_array_equal(c[key],c0[key]);np.testing.assert_array_equal(g[key],g0[key])
        for key in kw:
            if isinstance(kw[key],np.ndarray):np.testing.assert_array_equal(kw[key],before[key]);kw[key].fill(0)
        a['proposal']['velocity_m_s'][:]=999;a['correction']['temperature_K'][:]=999
        a['diagnostics']['omissions'].append('MUTATED');a['diagnostics']['free_axes'][:]=False
        again=p.apply(c,g)
        for key in FIELDS:np.testing.assert_array_equal(again['proposal'][key],saved['proposal'][key])
        self.assertNotIn('MUTATED',again['diagnostics']['omissions'])
        for name in ('matrix','mass_dof','Jv','Jr','JT','E','kappa','ids','free'):
            self.assertFalse(getattr(p,name).flags.writeable,name)

    def test_07_zero_defect_does_not_redefine_event_mean(self):
        c=fields(1);c['mean_velocity_m_s'][0,0]=-1/18;c['velocity_m_s'][0,0]=-2/3
        got=MaterialStageResponse(**scalar()).apply(c,copy.deepcopy(c))
        for key in FIELDS:
            np.testing.assert_array_equal(got['proposal'][key],c[key]);np.testing.assert_array_equal(got['correction'][key],np.zeros_like(c[key]))
        self.assertNotEqual(got['proposal']['mean_velocity_m_s'][0,0],1/6)
        for name in ('accepted','physical_map','restart_allowed','full_tangent'):self.assertIs(got['diagnostics'][name],False)

    def test_08_no_density_temperature_or_free_velocity_clipping(self):
        c=fields(1);g=copy.deepcopy(c);g['mean_velocity_m_s'][0,0]=100;g['temperature_K'][0]=-20
        r=MaterialStageResponse(**scalar()).apply(c,g)
        self.assertLess(r['proposal']['density_kg_m3'][0],0);self.assertEqual(r['proposal']['temperature_K'][0],-20)
        self.assertGreater(r['proposal']['mean_velocity_m_s'][0,0],20)
        np.testing.assert_array_equal(c['density_kg_m3'],[5.])

    def test_09_invalid_and_original_continuity_domain_fail_closed(self):
        changes=[('mass',[-1.]),('volume',[0.]),('material_tangent',np.full((1,9,9),np.nan)),
            ('material_temperature',np.full((1,9),np.inf)),('gradient_map',np.zeros((1,8,3))),
            ('rho_reference',0.),('rho_parent',[-1.]),('dt',0.),('bulk',0.),
            ('divergence',[-1.]),('divergence',[1.]),('free_axes',[[.5,1,1]]),('particle_ids',[1.5])]
        for key,value in changes:
            with self.subTest(key=key,value=str(value)[:20]):
                kw=scalar();kw[key]=value
                with self.assertRaises((ValueError,TypeError)):MaterialStageResponse(**kw)
        kw=coupled();kw['particle_ids'][:]=17
        with self.assertRaises(ValueError):MaterialStageResponse(**kw)
        p=MaterialStageResponse(**scalar());c=fields(1)
        for key in FIELDS:
            g=copy.deepcopy(c);g[key].flat[0]=np.nan
            with self.assertRaises(ValueError):p.apply(c,g)

    def test_10_exact_singularity_rejects_without_pinv(self):
        kw=scalar();kw.update(mass=np.ones(1),volume=np.ones(1),dt=1.,bulk=1.,rho_reference=1.,rho_parent=np.ones(1))
        kw['gradient_map'][:]=0;kw['gradient_map'][0,1,0]=1
        kw['material_tangent'][0,1,1]=-4
        # Pure Lxy has zero divergence. Literal Axx=1+.25*(-4)=0.
        p=MaterialStageResponse(**kw);c=fields(1);g=copy.deepcopy(c);g['mean_velocity_m_s'][0,0]=1
        with mock.patch('numpy.linalg.pinv',side_effect=AssertionError('pinv forbidden')):
            with self.assertRaises((ValueError,np.linalg.LinAlgError)):p.apply(c,g)

    def test_11_particle_permutation_and_id_values_do_not_change_action(self):
        kw=coupled();c,g=coupled_fields();base=MaterialStageResponse(**kw).apply(c,g)
        ids_only=copy.deepcopy(kw);ids_only['particle_ids']=np.array([900,3]);same=MaterialStageResponse(**ids_only).apply(c,g)
        for key in FIELDS:np.testing.assert_array_equal(base['correction'][key],same['correction'][key])
        order=np.array([1,0]);dofs=np.array([3,4,5,0,1,2]);q=copy.deepcopy(kw)
        for name in ('mass','volume','material_tangent','material_temperature','rho_parent','divergence','free_axes','particle_ids'):q[name]=q[name][order]
        q['gradient_map']=q['gradient_map'][order][:,:,dofs]
        rearranged=MaterialStageResponse(**q).apply({k:v[order] for k,v in c.items()},{k:v[order] for k,v in g.items()})
        for key in FIELDS:self.close(rearranged['correction'][key],base['correction'][key][order])

    def test_12_dimensional_unit_covariance_not_physical_scaling(self):
        # Manufactured coordinate-unit conversion only: L'=5L,t'=3t,m'=7m,T'=2T.
        kw=coupled();c,g=coupled_fields();base=MaterialStageResponse(**kw).apply(c,g);q=copy.deepcopy(kw)
        q['mass']*=7;q['volume']*=125;q['gradient_map']/=5;q['material_tangent']*=7/15
        q['material_temperature']*=7/90;q['bulk']*=7/45;q['rho_reference']*=7/125
        q['rho_parent']*=7/125;q['dt']*=3;q['divergence']/=3
        scales=dict(position_m=5.,velocity_m_s=5/3,mean_velocity_m_s=5/3,density_kg_m3=7/125,temperature_K=2.)
        result=MaterialStageResponse(**q).apply({k:v*scales[k] for k,v in c.items()},{k:v*scales[k] for k,v in g.items()})
        for key in FIELDS:self.close(result['correction'][key],base['correction'][key]*scales[key])

if __name__=='__main__':unittest.main(verbosity=2)
