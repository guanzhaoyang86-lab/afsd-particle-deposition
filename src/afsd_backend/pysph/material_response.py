"""Fixed-branch numerical response of the existing FKS/J2 material.

This is a search-direction component, not a material update or heat source.
The physical material.py, its 48 bisections, parameters and histories remain
unchanged. Yield/rate switching surfaces are explicitly outside this tangent.
"""
import numpy as np
from .material import FKS, rotation

CLASSIFICATION = 'FIXED_BRANCH_FKS_NUMERICAL_RESPONSE_NOT_PHYSICAL_UPDATE'


class NonDifferentiableMaterial(ValueError):
    """The requested classical material derivative is not defined/certified."""


def _real(value, name):
    a = np.asarray(value)
    if a.dtype.kind not in 'iuf' or not np.isfinite(a).all():
        raise ValueError('finite real '+name+' required')
    return np.array(a, dtype=float, copy=True)


def _tensor(value, name):
    a = _real(value, name)
    if a.ndim < 2 or a.shape[-2:] != (3, 3):
        raise ValueError(name+' must end in 3x3')
    return a


def _deviatoric(value, name):
    a = _tensor(value, name)
    scale = np.maximum(1., np.max(np.abs(a), axis=(-2, -1)))
    tolerance = 128*np.finfo(float).eps*scale
    if (np.any(np.max(np.abs(a-a.swapaxes(-2, -1)), axis=(-2, -1)) > tolerance)
            or np.any(np.abs(np.trace(a, axis1=-2, axis2=-1)) > 3*tolerance)):
        raise ValueError(name+' must be symmetric deviatoric; not projected')
    return a


def fks_partials(law, ep, rate, temperature):
    """Y, dY/dep [Pa], dY/dr [Pa s], dY/dT [Pa/K] on one FKS branch.

    ep is the END plastic strain of the already solved base point.
    The existing rate floor is represented, not moved or smoothed.
    """
    if not isinstance(law, FKS):
        raise ValueError('explicit existing FKS law required')
    attributes = ('a1','a2','a3','b1','b2','c1','c2','room','melt','reference_rate','rate_floor')
    values = np.array([getattr(law, k) for k in attributes], float)
    if (not np.isfinite(values).all() or law.a1 <= 0 or law.a2 < 0 or law.a3 <= 0
            or law.b1 < 0 or law.b2 < 1 or law.c1 <= 0 or law.c2 <= 0
            or law.melt <= law.room or law.reference_rate <= 0 or law.rate_floor <= 0):
        raise ValueError('unsupported FKS parameter domain')
    ep, r, T = np.broadcast_arrays(_real(ep,'ep'), _real(rate,'rate'), _real(temperature,'temperature'))
    if np.any(ep < 0) or np.any(r < 0) or np.any(T < law.room) or np.any(T >= law.melt):
        raise ValueError('FKS input outside unchanged solid domain')
    if np.any(r == law.rate_floor):
        raise NonDifferentiableMaterial('rate equals floor; no two-sided FKS derivative')
    span = law.melt-law.room
    ts = (T-law.room)/span
    above = r > law.rate_floor
    effective = np.where(above, r, law.rate_floor)
    log_rate = np.log(effective/law.reference_rate)
    lam = 1+law.b1*ts**law.b2*log_rate
    log_b = -np.logaddexp(0., -law.c1*ts)/law.c2
    theta = -np.expm1(log_b)
    hard = law.a1+law.a2*np.arctan(law.a3*ep)
    value = hard*lam*theta
    d_ep = law.a2*law.a3/(1+(law.a3*ep)**2)*lam*theta
    # effective is positive even below the floor; there is no invalid 0/0.
    d_rate = np.where(above, hard*theta*law.b1*ts**law.b2/effective, 0.)
    d_lam_T = law.b1*law.b2*ts**(law.b2-1)*log_rate/span
    d_theta_T = -(law.c1/law.c2)*np.exp(log_b)/(1+np.exp(law.c1*ts))/span
    d_T = hard*(theta*d_lam_T+lam*d_theta_T)
    if np.any(value <= 0) or not all(np.isfinite(x).all() for x in (value,d_ep,d_rate,d_T)):
        raise ValueError('nonpositive/nonfinite FKS response; no clipping')
    return dict(yield_stress_Pa=value.copy(), d_ep_Pa=d_ep.copy(),
                d_rate_Pa_s=d_rate.copy(), d_temperature_Pa_K=d_T.copy(),
                rate_branch=np.where(above,'ABOVE_FLOOR','BELOW_FLOOR'),
                classification=CLASSIFICATION)


def radial_direction(trial_stress, plastic_increment, delta_trial_stress,
                     delta_rate, delta_temperature, shear, partials, active):
    """Directional derivative of radial return at an already solved base.

    Not a fresh return solve. The caller must supply that base's increment and
    end-state Y partials. Base closure error is exposed, never turned into heat.
    Inactive qtr==0 and yield equality are rejected as norm/switching boundaries.
    """
    s = _deviatoric(trial_stress,'trial_stress')
    ds = _deviatoric(delta_trial_stress,'delta_trial_stress')
    dep = _real(plastic_increment,'plastic_increment')
    dr = _real(delta_rate,'delta_rate'); dT = _real(delta_temperature,'delta_temperature')
    G = _real(shear,'shear')
    mask = np.asarray(active)
    if mask.dtype.kind != 'b':
        raise ValueError('explicit boolean active branch required')
    keys = ('yield_stress_Pa','d_ep_Pa','d_rate_Pa_s','d_temperature_Pa_K')
    if not isinstance(partials,dict) or any(k not in partials for k in keys):
        raise ValueError('complete material partials required')
    y, ye, yr, yt = [_real(partials[k],k) for k in keys]
    shape = np.broadcast_shapes(s.shape[:-2],ds.shape[:-2],dep.shape,dr.shape,dT.shape,
                                G.shape,mask.shape,y.shape,ye.shape,yr.shape,yt.shape)
    s = np.broadcast_to(s,shape+(3,3)); ds = np.broadcast_to(ds,shape+(3,3))
    dep,dr,dT,G,mask,y,ye,yr,yt = [np.broadcast_to(a,shape) for a in (dep,dr,dT,G,mask,y,ye,yr,yt)]
    if np.any(G <= 0) or np.any(dep < 0) or np.any(y <= 0):
        raise ValueError('positive modulus/yield and nonnegative base increment required')
    if np.any(mask & (dep <= 0)) or np.any((~mask) & (dep != 0)):
        raise ValueError('base increment inconsistent with explicit branch')
    qtr = np.sqrt(1.5*np.sum(s*s,axis=(-2,-1)))
    if np.any(qtr == 0):
        raise NonDifferentiableMaterial('equivalent-stress norm at zero is not differentiable')
    if np.any((~mask) & (qtr == y)):
        raise NonDifferentiableMaterial('yield boundary is not a certified fixed branch')
    if np.any((~mask) & (qtr > y)):
        raise ValueError('inactive branch is outside yield')
    q = qtr-3*G*dep
    denominator = 3*G+ye
    if np.any(q <= 0) or np.any(denominator <= 0):
        raise ValueError('nonpositive returned stress or implicit return denominator')
    dqtr = 1.5*np.sum(s*ds,axis=(-2,-1))/qtr
    ddep = np.where(mask,(dqtr-yr*dr-yt*dT)/denominator,0.)
    dq = dqtr-3*G*ddep
    scale = q/qtr
    dscale = dq/qtr-q*dqtr/(qtr*qtr)
    result = scale[...,None,None]*ds+dscale[...,None,None]*s
    if not all(np.isfinite(a).all() for a in (result,ddep,dq)):
        raise ValueError('nonfinite material response; no clipping')
    return dict(delta_stress_Pa=result.copy(),delta_plastic_strain=ddep.copy(),
                delta_equivalent_Pa=dq.copy(),base_return_residual_Pa=np.where(mask,q-y,0.),
                classification=CLASSIFICATION,physical_update=False)


def trial_direction(stress, gradient, delta_gradient, dt, shear):
    """Analytic Frechet derivative of the unchanged corotational trial stress.

    Rotations use the same Rodrigues baseline as material.update. Series only
    avoid cancellation in derivatives near zero spin; no physical field is
    projected or clipped. The returned rate derivative is undefined at zero
    deviatoric rate for nonzero deviatoric directions and explicitly rejects.
    """
    s = _deviatoric(stress,'stress')
    L = _tensor(gradient,'gradient'); dL = _tensor(delta_gradient,'delta_gradient')
    h = _real(dt,'dt'); G = _real(shear,'shear')
    if h.shape or G.shape or h <= 0 or G <= 0:
        raise ValueError('positive scalar dt/shear required')
    shape = np.broadcast_shapes(s.shape[:-2],L.shape[:-2],dL.shape[:-2])
    s,L,dL = [np.broadcast_to(a,shape+(3,3)) for a in (s,L,dL)]
    eye = np.eye(3)
    def dev(a):
        return .5*(a+a.swapaxes(-2,-1))-np.trace(a,axis1=-2,axis2=-1)[...,None,None]*eye/3
    D,dD = dev(L),dev(dL)
    A = .25*float(h)*(L-L.swapaxes(-2,-1))
    dA = .25*float(h)*(dL-dL.swapaxes(-2,-1))
    squared = .5*np.sum(A*A,axis=(-2,-1)); angle = np.sqrt(squared)
    d_squared = np.sum(A*dA,axis=(-2,-1))
    small = squared < 1e-6
    # Evaluate the non-small denominator at 1 for series nodes, avoiding 0/0.
    angle_safe = np.where(small,1.,angle)
    da_large = (angle_safe*np.cos(angle_safe)-np.sin(angle_safe))/(2*angle_safe**3)
    db_large = (angle_safe*np.sin(angle_safe)-4*np.sin(angle_safe/2)**2)/(2*angle_safe**4)
    da = np.where(small,-1/6+squared/60-squared**2/1680+squared**3/90720,da_large)
    db = np.where(small,-1/24+squared/360-squared**2/13440+squared**3/907200,db_large)
    a = np.sinc(angle/np.pi); b = .5*np.sinc(angle/(2*np.pi))**2
    R = rotation(A)
    dR = ((da*d_squared)[...,None,None]*A+a[...,None,None]*dA
          +(db*d_squared)[...,None,None]*(A@A)+b[...,None,None]*(dA@A+A@dA))
    Q = R@R; dQ = dR@R+R@dR
    Rt=R.swapaxes(-2,-1); dRt=dR.swapaxes(-2,-1)
    Qt=Q.swapaxes(-2,-1); dQt=dQ.swapaxes(-2,-1)
    trial = Q@s@Qt+2*float(G*h)*(R@D@Rt)
    delta = dQ@s@Qt+Q@s@dQt+2*float(G*h)*(dR@D@Rt+R@dD@Rt+R@D@dRt)
    rate = np.sqrt(2/3*np.sum(D*D,axis=(-2,-1)))
    zero = rate == 0
    if np.any(zero & np.any(dD != 0,axis=(-2,-1))):
        raise NonDifferentiableMaterial('rate norm at zero with nonzero deviatoric direction')
    dr = np.divide((2/3)*np.sum(D*dD,axis=(-2,-1)),rate,
                   out=np.zeros(shape),where=~zero)
    if not all(np.isfinite(a).all() for a in (trial,delta,rate,dr)):
        raise ValueError('nonfinite corotational response')
    return dict(trial_stress_Pa=trial,delta_trial_stress_Pa=delta,
                rate_per_s=rate,delta_rate_per_s=dr,classification=CLASSIFICATION)

