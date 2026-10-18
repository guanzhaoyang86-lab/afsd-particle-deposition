"""3D corotational hypoelastic J2/FKS update; no Johnson-Cook fallback."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class FKS:
    a1: float = 317.6e6
    a2: float = 1e6
    a3: float = 70.0
    b1: float = 12.44  # DEBUG dimensionless interpretation, NOT MPa conversion.
    b2: float = 5.65
    c1: float = 13.0
    c2: float = 0.01
    room: float = 293.15
    melt: float = 855.15
    reference_rate: float = 1.0
    rate_floor: float = 1.0  # Explicit DEBUG low-rate regularization.

    def yield_stress(self, ep, rate, temperature):
        ep, rate, temperature = np.broadcast_arrays(ep, rate, temperature)
        if not all(np.isfinite(a).all() for a in (ep, rate, temperature)):
            raise ValueError('non-finite FKS input')
        if np.any(ep < 0) or np.any(rate < 0):
            raise ValueError('negative strain or strain rate')
        if np.any(temperature < self.room) or np.any(temperature >= self.melt):
            raise ValueError('FKS temperature outside explicitly supported solid domain')
        ts = (temperature-self.room)/(self.melt-self.room)
        lam = 1+self.b1*ts**self.b2*np.log(np.maximum(rate, self.rate_floor)/self.reference_rate)
        # 1 - (1 + exp(-c1*T*))**(-1/c2), stable near high T.
        theta = -np.expm1(-np.logaddexp(0, -self.c1*ts)/self.c2)
        result = (self.a1+self.a2*np.arctan(self.a3*ep))*lam*theta
        if np.any(result <= 0) or not np.isfinite(result).all():
            raise ValueError('nonpositive/non-finite FKS yield; not clipped')
        return result


def equivalent(stress):
    return np.sqrt(1.5*np.sum(stress*stress, axis=(-2, -1)))


def rotation(spin_dt):
    """Batched Rodrigues exponential, exact for constant antisymmetric spin."""
    w = np.asarray(spin_dt)
    angle = np.sqrt(np.sum(w*w, axis=(-2, -1))/2)
    a = np.sinc(angle/np.pi)
    b = 0.5*np.sinc(angle/(2*np.pi))**2
    return np.eye(3)+a[..., None, None]*w+b[..., None, None]*(w@w)


def update(stress, ep, gradient, temperature, dt, shear, law):
    if dt <= 0 or shear <= 0:
        raise ValueError('positive dt and shear modulus required')
    d = 0.5*(gradient+gradient.swapaxes(-1, -2))
    dev = d-np.trace(d, axis1=-2, axis2=-1)[..., None, None]*np.eye(3)/3
    spin = 0.5*(gradient-gradient.swapaxes(-1, -2))
    half = rotation(spin*(dt/2))
    full = half@half
    trial = full@stress@full.swapaxes(-1, -2)+2*shear*dt*(half@dev@half.swapaxes(-1, -2))
    qtrial = equivalent(trial)
    rate = np.sqrt(2/3*np.sum(dev*dev, axis=(-2, -1)))
    active = qtrial > law.yield_stress(ep, rate, temperature)
    lo = np.zeros_like(ep)
    hi = np.where(active, qtrial/(3*shear), 0)
    for _ in range(48):
        mid = (lo+hi)/2
        residual = qtrial-3*shear*mid-law.yield_stress(ep+mid, rate, temperature)
        lo = np.where(active & (residual > 0), mid, lo)
        hi = np.where(active & (residual <= 0), mid, hi)
    dep = (lo+hi)/2
    q = qtrial-3*shear*dep
    scale = np.divide(q, qtrial, out=np.ones_like(q), where=qtrial > 0)
    out = trial*scale[..., None, None]
    # Backward-Euler plastic work sigma_new:deps_p = q_new * dep.
    # Algorithmic return dissipation 1.5 G dep^2 is reported separately.
    return out, ep+dep, q*dep, 1.5*shear*dep*dep, rate
