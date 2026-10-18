"""P0 page 5 Eqs 2-3, visually checked; three-dimensional only."""
import numpy as np


class HyperbolicSpline:
    dim = 3
    radius_scale = 2.0

    def value(self, r, h):
        q = np.asarray(r) / h
        if np.any(q < 0) or np.any(np.asarray(h) <= 0):
            raise ValueError('negative radius or nonpositive smoothing length')
        f = np.where(q < 1, q**3 - 6*q + 6, np.where(q < 2, (2-q)**3, 0))
        return 15 / (62*np.pi*h**3) * f

    def gradient(self, displacement, h):
        x = np.asarray(displacement, dtype=float)
        r = np.linalg.norm(x, axis=-1)
        if np.any(r == 0):
            raise ValueError('distinct coincident particles: cusp gradient undefined')
        q = r/h
        df = np.where(q < 1, 3*q*q-6, np.where(q < 2, -3*(2-q)**2, 0))
        return (15/(62*np.pi*h**4)*df/r)[..., None]*x

    def self_gradient(self):
        # Isotropic self contribution is excluded, not a differentiable origin.
        return np.zeros(3)
