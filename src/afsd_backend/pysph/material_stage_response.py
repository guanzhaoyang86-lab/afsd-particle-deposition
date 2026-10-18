"""Journal-frozen rate-aware material search block, never a physical map.

It retains the existing five-field correction interface. Geometry, contact
response and BE source feedback are NOT differentiated. The actual residual
must still be evaluated by the unchanged process trial.
"""
import copy,hashlib,json
import numpy as np
from .material import FKS
from .material_response import fks_partials,radial_direction,trial_direction,NonDifferentiableMaterial
from .stage_preconditioner import _fields,MAXIMUM_BACKWARD_ERROR
from .stage_iteration import FIELDS

CLASSIFICATION='JOURNAL_FROZEN_RATE_AWARE_MATERIAL_SEARCH_NOT_PHYSICAL_MAP'


def _array(value,shape,name):
    a=np.asarray(value)
    if a.dtype.kind not in 'iuf' or a.shape!=shape or not np.isfinite(a).all():
        raise ValueError('invalid '+name)
    return np.array(a,float,copy=True)


def _close(a,b,atol,rtol,name):
    a,b=np.asarray(a),np.asarray(b)
    if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(abs(a-b)>atol+rtol*abs(b)):
        raise ValueError('journal mismatch '+name)


def _fingerprint(value):
    def plain(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,dict):return {str(k):plain(q) for k,q in v.items()}
        if isinstance(v,(list,tuple)):return [plain(q) for q in v]
        return v
    return hashlib.sha256(json.dumps(plain(value),sort_keys=True,allow_nan=False).encode()).hexdigest()


class MaterialStageResponse:
    """Eliminate density in a frozen local material approximation.

    C_alg is dS1/dL [Pa s], not a modulus. Returned reactions are only algebraic.
    The full generally nonsymmetric A is kept; no symmetrization, pinv or mass
    modification. Full constrained columns participate in pressure backsolve.
    """
    def __init__(self,mass,volume,gradient_map,material_tangent,material_temperature,
                 bulk,rho_reference,rho_parent,dt,divergence,free_axes,particle_ids):
        ids=np.asarray(particle_ids);n=len(ids)
        if ids.shape!=(n,) or ids.dtype.kind not in 'iu' or n<1 or len(np.unique(ids))!=n:
            raise ValueError('unique integer particle IDs')
        self.ids=ids.copy();self.count=n
        m=_array(mass,(n,),'mass');V=_array(volume,(n,),'volume')
        D=_array(gradient_map,(n,9,3*n),'gradient_map')
        C=_array(material_tangent,(n,9,9),'material_tangent')
        ST=_array(material_temperature,(n,9),'material_temperature')
        rhoref=_array(np.broadcast_to(rho_reference,(n,)),(n,),'reference_density')
        rhoold=_array(rho_parent,(n,),'parent_density')
        div=_array(divergence,(n,),'divergence')
        h=_array(dt,(),'dt');K=_array(bulk,(),'bulk')
        free=np.asarray(free_axes)
        if free.shape!=(n,3) or free.dtype.kind not in 'biuf' or not np.isfinite(free).all() or not np.all((free==0)|(free==1)):
            raise ValueError('exact free-axis mask')
        if np.any(m<=0) or np.any(V<=0) or np.any(rhoref<=0) or np.any(rhoold<=0) or h<=0 or K<=0:
            raise ValueError('positive physical mass,volume,density,dt,bulk')
        denominator=1+.5*h*div
        if np.any(denominator<=0) or np.any(1-.5*h*div<=0):
            raise ValueError('original continuity domain')
        E=np.trace(D.reshape(n,3,3,3*n),axis1=1,axis2=2)
        Jv=-.5*np.einsum('nik,nij,njl->kl',D,V[:,None,None]*C,D,optimize=True)
        Jr=.5*E.T*(V*float(K)/rhoref)[None,:]
        JT=-.25*np.einsum('nik,ni->kn',D,V[:,None]*ST,optimize=True)
        kappa=-float(h)*rhoold/denominator**2
        mass_dof=np.repeat(m,3)
        A=np.diag(mass_dof)-.5*float(h)*Jv-.5*float(h)*((Jr*kappa[None,:])@E)
        if not np.isfinite(A).all():
            raise ValueError('nonfinite material search matrix')
        self.dt=float(h);self.free=free.astype(bool).ravel();self.matrix=A
        self.mass_dof=mass_dof;self.Jr=Jr;self.JT=JT;self.kappa=kappa;self.E=E
        self.Jv=Jv
        self.metadata=dict(classification=CLASSIFICATION,accepted=False,physical_map=False,
            restart_allowed=False,full_tangent=False,dt_s=self.dt,
            omissions=['geometry_derivatives','contact_response','BE_heat_feedback'],
            continuity='EXACT_CANDIDATE_DENSITY_SCHUR_IN_NUMERICAL_MODEL',
            temperature='DIRECT_DEFECT_WITH_MATERIAL_TMID_FORCE_COUPLING_ONLY')
        for a in (self.ids,self.free,self.matrix,self.mass_dof,self.Jr,self.JT,self.kappa,self.E,self.Jv):
            a.flags.writeable=False

    @classmethod
    def from_journal(cls,parent,state,p,candidate,evaluation):
        """Assemble only from an explicitly bound saved or current trial journal."""
        e=evaluation;ids=np.asarray(parent['particle_id']);n=len(ids)
        c=_fields(candidate,n)
        if e.get('accepted') is not False or e.get('restart_allowed') is not False:
            raise ValueError('uncommitted evaluation required')
        if not np.array_equal(ids,np.asarray(e['particle_id'])) or e['time_before_s']!=state['t']:
            raise ValueError('journal parent time/ID mismatch')
        if _fingerprint(state['parameters'])!=_fingerprint(p):
            raise ValueError('journal physical parameters differ from parent')
        h=float(e['dt_s']);G=float(p['shear']);K=2*G*(1+p['poisson'])/(3*(1-2*p['poisson']))
        if h<=0 or not np.isfinite([h,G,K]).all():
            raise ValueError('positive finite interval/moduli required')
        mid=e['midpoint'];mat=e['material']
        mass=_array(parent['m'],(n,),'parent_mass');rho=_array(parent['rho'],(n,),'parent_density')
        if np.any(mass<=0) or np.any(rho<=0):raise ValueError('positive parent mass/density')
        ep0=_array(parent['ep'],(n,),'parent_ep');S0=np.asarray(parent['S'],float).reshape(n,3,3)
        oldx=np.column_stack([parent[k] for k in ('x','y','z')])
        _close(mid['mass_kg'],mass,0,0,'mass')
        _close(mid['position_m'],.5*(oldx+c['position_m']),0,0,'midpoint_position')
        _close(mid['density_kg_m3'],.5*(rho+c['density_kg_m3']),0,0,'midpoint_density')
        _close(mid['temperature_K'],.5*(np.asarray(parent['temperature'])+c['temperature_K']),0,0,'midpoint_temperature')
        _close(mid['mean_velocity_m_s'],c['mean_velocity_m_s'],0,0,'independent_mean_velocity')
        for key in FIELDS:
            _close(e['residual'][key],c[key]-np.asarray(e['map'][key]),0,0,'raw_residual_'+key)
        rhomid=np.asarray(mid['density_kg_m3'],float)
        if np.any(rhomid<=0):raise ValueError('positive midpoint density')
        V=mass/rhomid
        i=np.asarray(mid['i']);j=np.asarray(mid['j'])
        if i.ndim!=1 or j.shape!=i.shape or i.dtype.kind not in 'iu' or j.dtype.kind not in 'iu' or np.any(i<0) or np.any(j<0) or np.any(i>=n) or np.any(j>=n) or np.any(i==j):
            raise ValueError('invalid saved pair identities')
        B=_array(mid['B'],(n,3,3),'B');g=_array(mid['g'],(len(i),3),'g')
        gi=np.einsum('nji,nj->ni',B[i],g);gj=np.einsum('nji,nj->ni',B[j],g)
        D=np.zeros((n,3,3,3*n))
        for component in range(3):
            for axis in range(3):
                np.add.at(D,(i,component,axis,3*j+component),V[j]*gi[:,axis])
                np.add.at(D,(i,component,axis,3*i+component),-V[j]*gi[:,axis])
                np.add.at(D,(j,component,axis,3*j+component),V[i]*gj[:,axis])
                np.add.at(D,(j,component,axis,3*i+component),-V[i]*gj[:,axis])
        D=D.reshape(n,9,3*n);L=_array(mid['velocity_gradient_per_s'],(n,3,3),'L')
        rebuilt=np.einsum('nij,j->ni',D,c['mean_velocity_m_s'].ravel()).reshape(n,3,3)
        _close(rebuilt,L,1e-10,1e-12,'saved_gradient_operator_identity')
        law=FKS(**{k:p[k] for k in FKS.__dataclass_fields__})
        ep1=_array(mat['ep'],(n,),'returned_ep');dep=ep1-ep0
        if np.any(dep<0):raise ValueError('decreasing plastic history')
        active=dep>0
        base=trial_direction(S0,L,np.zeros_like(L),h,G)
        rate=base['rate_per_s'];trial=base['trial_stress_Pa']
        _close(mat['strain_rate_per_s'],rate,1e-12,1e-12,'strain_rate')
        partials=fks_partials(law,ep1,rate,mid['temperature_K'])
        initial=fks_partials(law,ep0,rate,mid['temperature_K'])
        qtr=np.sqrt(1.5*np.sum(trial*trial,axis=(1,2)))
        if np.any(qtr==0):raise NonDifferentiableMaterial('zero equivalent trial stress')
        if np.any(qtr==initial['yield_stress_Pa']):
            raise NonDifferentiableMaterial('journal at yield switching surface')
        if not np.array_equal(active,qtr>initial['yield_stress_Pa']):
            raise ValueError('saved active branch inconsistent')
        returned=(1-3*G*dep/qtr)[:,None,None]*trial
        _close(mat['deviatoric_stress_Pa'],returned,1e-6,1e-12,'saved_returned_stress')
        C=np.zeros((n,9,9))
        for a in range(3):
            for b in range(3):
                basis=np.zeros((n,3,3));basis[:,a,b]=1.
                tr=trial_direction(S0,L,basis,h,G)
                answer=radial_direction(trial,dep,tr['delta_trial_stress_Pa'],
                    tr['delta_rate_per_s'],np.zeros(n),G,partials,active)
                C[:,:,3*a+b]=answer['delta_stress_Pa'].reshape(n,9)
        temp=radial_direction(trial,dep,np.zeros_like(trial),np.zeros(n),np.ones(n),G,partials,active)
        traces=e['boundary']['traces']
        if len(traces)!=n or any(tr['particle_index']!=k or tr['particle_id']!=ids[k] for k,tr in enumerate(traces)):
            raise ValueError('constraint trace IDs/order')
        free=np.asarray([tr['free_axes'] for tr in traces])
        obj=cls(mass,V,D,C,temp['delta_stress_Pa'].reshape(n,9),K,p['rho0'],rho,h,
                np.trace(L,axis1=1,axis2=2),free,ids)
        obj.metadata.update(parent_time_s=state['t'],parent_evaluation_hash=e['accepted_state_hash'],
            base_candidate_sha256=_fingerprint(c),base_evaluation_sha256=_fingerprint(e),
            plastic_particle_ids=ids[active].copy(),pair_count=len(i),
            material_base_return_residual_Pa=temp['base_return_residual_Pa'].copy())
        return obj

    def apply(self,candidate,mapping):
        c=_fields(candidate,self.count);mapped=_fields(mapping,self.count)
        defect={k:mapped[k]-c[k] for k in FIELDS}
        dbar=defect['mean_velocity_m_s'].ravel()
        rhs=self.mass_dof*dbar+.5*self.dt*(self.Jr@defect['density_kg_m3']+self.JT@defect['temperature_K'])
        dv=np.zeros_like(dbar);free=self.free;dv[~free]=dbar[~free]
        A=self.matrix;aff=A[np.ix_(free,free)]
        shifted=rhs[free]-A[np.ix_(free,~free)]@dv[~free]
        if np.any(free):dv[free]=np.linalg.solve(aff,shifted)
        residual=aff@dv[free]-shifted
        norm=lambda x:float(np.max(abs(x),initial=0.))
        matrix_norm=float(np.max(np.sum(abs(aff),axis=1),initial=0.))
        denominator=matrix_norm*norm(dv[free])+norm(shifted)
        backward=norm(residual)/denominator if denominator else 0.
        if not np.isfinite(dv).all() or not np.isfinite(backward) or backward>MAXIMUM_BACKWARD_ERROR:
            raise ValueError('material numerical solve failed original1e-12 backward error')
        change=(dv-dbar).reshape(self.count,3)
        correction={k:v.copy() for k,v in defect.items()}
        correction['mean_velocity_m_s']=dv.reshape(self.count,3).copy()
        correction['position_m']+=self.dt*change
        correction['velocity_m_s']+=2*change
        correction['density_kg_m3']+=self.kappa*(self.E@dv)
        proposal={k:c[k]+correction[k] for k in FIELDS}
        if not all(np.isfinite(a).all() for a in proposal.values()):
            raise ValueError('nonfinite numerical proposal; no clipping')
        return dict(proposal=proposal,correction=correction,diagnostics=dict(copy.deepcopy(self.metadata),
            free_axes=free.reshape(self.count,3).copy(),backward_error_inf=backward,
            free_linear_residual_Ns=residual.copy(),algebraic_reaction_impulse_Ns=A@dv-rhs,
            reaction_scope='SEARCH_ALGEBRA_NOT_PHYSICAL_REACTION',
            matrix_asymmetry_inf_kg=norm(A-A.T)))
