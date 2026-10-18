"""NumPy batching of the frozen finite-patch geometry, with unchanged contact.

Only endpoint/midpoint geometry is batched. The scalar ``_patches`` and
``_sweep`` are inherited unchanged, including actual swept-solid rejection.
DG forces, spring return, signed work, and every acceptance check are inherited
from ChainRuleDiscreteContact. This is an implementation variant, not a new
physical model, event integrator, broad-phase approximation, or history cache.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))

from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from afsd_backend.pysph.quasistatic_contact import _pose, ContactGeometryError


IMPLEMENTATION = 'NUMPY_BATCHED_FINITE_PATCH_EQUIVALENT'


class BatchedChainRuleContact(ChainRuleDiscreteContact):
    implementation = IMPLEMENTATION

    def _patches_batch(self, points, reject_active_ties=True):
        """Same seven projections and local tie bound as scalar ``_patches``.

        Both annular caps retain their separate squared-distance argmin; no
        nearest-face pruning or averaging across a medial-axis tie is used.
        Rows within an arithmetic bound of a geometric/contact branch, as
        well as every invalid row, invoke the inherited scalar oracle in row
        order. This preserves scalar accept/reject decisions even where a
        batched norm differs in its last bit from a one-vector norm.
        """
        points = np.asarray(points, float)
        a, R, H = self.half_bore, self.radius, self.height
        n = len(points)
        xy = points[:, :2]
        radial = np.linalg.norm(xy, axis=1)
        largest_xy = np.max(np.abs(xy), axis=1)
        inside = ((points[:, 2] > 0.) & (points[:, 2] < H) &
                  (radial < R) & (largest_xy > a))
        circle = np.empty_like(xy)
        nonzero = radial > 0.
        circle[nonzero] = xy[nonzero] * (R / radial[nonzero, None])
        circle[~nonzero] = [R, 0.]

        closest = np.empty((n, 7, 3))
        ambiguity = np.zeros((n, 7), bool)
        outside_circle = radial > R
        inside_square = ~outside_circle & (largest_xy < a)
        # Outside square / inside circle: the planar point itself is nearest.
        planar = np.where(outside_circle[:, None], circle, xy)
        for feature, z in enumerate((0., H)):
            closest[:, feature, :2] = planar
            closest[:, feature, 2] = z

        rows = np.flatnonzero(inside_square)
        if len(rows):
            p = points[rows]
            options = np.repeat(p[:, None, :], 4, axis=1)
            # Exact scalar candidate order: X-, X+, Y-, Y+.
            options[:, 0, 0] = -a
            options[:, 1, 0] = a
            options[:, 2, 1] = -a
            options[:, 3, 1] = a
            for feature, z in enumerate((0., H)):
                options[:, :, 2] = z
                delta = p[:, None, :] - options
                d2 = np.einsum('nki,nki->nk', delta, delta)
                best = np.argmin(d2, axis=1)
                selected = options[np.arange(len(rows)), best]
                # Same cancellation certificate as scalar implementation;
                # adding a global distance tolerance would alter tie branches.
                terms = (options - selected[:, None, :]) * (
                    (options - p[:, None, :]) +
                    (selected - p)[:, None, :])
                tied = (np.abs(np.sum(terms, axis=2)) <=
                        64 * np.finfo(float).eps * np.sum(np.abs(terms), axis=2))
                distinct = np.any(options != selected[:, None, :], axis=2)
                closest[rows, feature] = selected
                ambiguity[rows, feature] = np.any(tied & distinct, axis=1)

        for axis in (0, 1):
            for side, sign in enumerate((-1., 1.)):
                feature = 2 + 2 * axis + side
                closest[:, feature] = points
                closest[:, feature, axis] = sign * a
                closest[:, feature, 1-axis] = np.clip(points[:, 1-axis], -a, a)
                closest[:, feature, 2] = np.clip(points[:, 2], 0., H)
        closest[:, 6, :2] = circle
        closest[:, 6, 2] = np.clip(points[:, 2], 0., H)
        ambiguity[:, 6] = radial == 0.
        delta = points[:, None, :] - closest
        distances = np.linalg.norm(delta, axis=2)
        invalid = inside | np.any(distances == 0., axis=1)
        if reject_active_ties:
            invalid |= np.any((distances < self.particle_radius) & ambiguity, axis=1)
        eps_bound = 128*np.finfo(float).eps
        length_scale = np.maximum(np.max(np.abs(points), axis=1), max(a, R, H))
        near_branch = ((np.abs(radial-R) <= eps_bound*length_scale) |
                       (np.abs(largest_xy-a) <= eps_bound*length_scale) |
                       (np.abs(points[:, 2]) <= eps_bound*length_scale) |
                       (np.abs(points[:, 2]-H) <= eps_bound*length_scale) |
                       np.any(distances <= eps_bound*length_scale[:, None], axis=1) |
                       np.any(np.abs(distances-self.particle_radius) <=
                              eps_bound*length_scale[:, None], axis=1))
        normals = np.zeros_like(delta)
        np.divide(delta, distances[:, :, None], out=normals,
                  where=distances[:, :, None] != 0.)
        for row in np.flatnonzero(invalid | near_branch):
            distances[row], normals[row], closest[row] = self._patches(
                points[row], reject_active_ties=reject_active_ties)
        return distances, normals, closest

    def geometry(self, positions, pose):
        x = np.asarray(positions, float)
        if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
            raise ValueError('finite (N,3) positions required')
        origin, rotation = _pose(pose)
        local = (x-origin) @ rotation
        gaps = np.full((len(x), 8), np.inf)
        normal_body = np.zeros((len(x), 8, 3))
        closest_body = np.zeros_like(normal_body)
        gaps[:, :7], normal_body[:, :7], closest_body[:, :7] = self._patches_batch(local)
        if self.base_z is not None:
            gaps[:, 7] = x[:, 2]-self.base_z
            if np.any(gaps[:, 7] < 0.):
                raise ContactGeometryError('particle centre below rigid base')
            normal_body[:, 7, 2] = 1.
            closest_body[:, 7] = x
            closest_body[:, 7, 2] = self.base_z
        active = gaps < self.particle_radius
        overlap = np.where(active, self.particle_radius-gaps, 0.)
        magnitude = self.kn*overlap
        normal_world = normal_body.copy()
        normal_world[:, :7] = normal_body[:, :7] @ rotation.T
        closest_world = closest_body.copy()
        closest_world[:, :7] = closest_body[:, :7] @ rotation.T+origin
        return dict(local_positions_m=local, gaps_m=gaps, active=active,
                    normal_body=normal_body, normal_world=normal_world,
                    closest_points_m=closest_world, magnitude_N=magnitude,
                    force_by_feature_N=magnitude[:, :, None]*normal_world,
                    potential_J=.5*self.kn*overlap**2, overlap_m=overlap)

    def _midpoint_geometry(self, local0, local1, x0, x1):
        middle = .5*(local0+local1)
        gaps = np.full((len(middle), 8), np.inf)
        normals = np.zeros((len(middle), 8, 3))
        gaps[:, :7], normals[:, :7], _ = self._patches_batch(middle)
        if self.base_z is not None:
            gaps[:, 7] = .5*(x0[:, 2]+x1[:, 2])-self.base_z
            normals[:, 7, 2] = 1.
        return gaps, normals
