"""Independent finite-patch, incremental-spring contact approximation.

This is not the inertial contact kick. Each of seven finite tool patches and
the base has its OWN distance penalty. In particular adjacent patches both
contribute at a lip. The resulting potential is conservative but is NOT the
old one-body nearest-distance potential. Particle centres must remain outside
the tool and above the base throughout the declared interpolation path.

Tangential history is a displacement in tool body coordinates (base: laboratory
coordinates). Coulomb radial return eliminates the tangential force unknown.
History projection/release and backward-Euler dissipation are reported apart
from physical sliding heat. No history or input array is mutated.
"""
import numpy as np


CLASSIFICATION = 'HOT_START_QUASISTATIC_FINITE_PATCH_CONTACT_APPROXIMATE_COMPONENT'
FEATURES = ('SHOULDER', 'TOP', 'BORE_X_MINUS', 'BORE_X_PLUS',
            'BORE_Y_MINUS', 'BORE_Y_PLUS', 'OUTER_CYLINDER', 'BASE')


class ContactGeometryError(ValueError):
    """Candidate rejected; a closest-point branch or admissibility is invalid."""


def _pose(pose):
    origin = np.asarray(pose['origin'], float)
    rotation = np.asarray(pose['rotation'], float)
    if (origin.shape != (3,) or rotation.shape != (3, 3) or
            not np.isfinite(origin).all() or not np.isfinite(rotation).all() or
            not np.allclose(rotation.T@rotation, np.eye(3), atol=1e-12, rtol=0.) or
            abs(np.linalg.det(rotation)-1.) > 1e-12):
        raise ValueError('finite origin and proper orthonormal rotation required')
    return origin, rotation


def _rotation_path(rotation0, rotation1):
    """SO(3) shortest geodesic; exactly pi is rejected as ambiguous."""
    relative = rotation0.T@rotation1
    skew = .5*np.array([relative[2, 1]-relative[1, 2],
                       relative[0, 2]-relative[2, 0],
                       relative[1, 0]-relative[0, 1]])
    sine = float(np.linalg.norm(skew))
    cosine = .5*(float(np.trace(relative))-1.)
    angle = float(np.arctan2(sine, cosine))
    if np.pi-angle <= 1e-12:
        raise ContactGeometryError('ambiguous pi tool rotation path')
    if angle == 0.:
        return angle, lambda fraction: rotation0
    axis = skew/sine
    K = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]],
                  [-axis[1], axis[0], 0.]])
    def rotation_at(fraction):
        theta = angle*fraction
        return rotation0@(np.eye(3)+np.sin(theta)*K+(1.-np.cos(theta))*(K@K))
    return angle, rotation_at


class QuasistaticContact:
    """SI inputs; stiffnesses N/m, dimensions m, coefficient dimensionless.

    ``evaluate`` accepts x0/x1 (N,3), poses with origin/rotation, dt in seconds,
    and history (N,8,3) in metres. Returned ``history`` is a candidate only.
    Inputs have no particle identity, tracer, AS/RS, mass, or temperature.
    ``base_z=None`` explicitly disables the base for isolated component tests.
    """
    features = FEATURES

    def __init__(self, tool_radius, bore_width, tool_height, particle_radius,
                 normal_stiffness, tangential_stiffness, friction, base_z=0.):
        values = np.asarray([tool_radius, bore_width, tool_height, particle_radius,
                             normal_stiffness, tangential_stiffness, friction], float)
        if (not np.isfinite(values).all() or np.any(values[:6] <= 0.) or
                friction < 0. or tool_radius <= bore_width/np.sqrt(2.) or
                particle_radius >= bore_width/2.):
            raise ValueError('positive geometry/stiffness, nonnegative friction required')
        if base_z is not None and not np.isfinite(base_z):
            raise ValueError('finite base_z or None required')
        self.radius = float(tool_radius)
        self.half_bore = .5*float(bore_width)
        self.height = float(tool_height)
        self.particle_radius = float(particle_radius)
        self.kn = float(normal_stiffness)
        self.kt = float(tangential_stiffness)
        self.mu = float(friction)
        self.base_z = None if base_z is None else float(base_z)

    def _patches(self, point, reject_active_ties=True):
        """Exact distance/projection to finite annular faces/rectangles/cylinder.

        Planar annular faces are a disk with an open square removed. At an
        active medial-axis tie the gradient does not exist, so reject instead
        of averaging normals. Inactive ties cannot contribute any force.
        """
        a, R, H = self.half_bore, self.radius, self.height
        radial = float(np.linalg.norm(point[:2]))
        if (0. < point[2] < H and radial < R and
                np.max(np.abs(point[:2])) > a):
            raise ContactGeometryError('particle centre inside finite tool solid')
        circle = point[:2]*(R/radial) if radial > 0. else np.array([R, 0.])
        if radial > R:
            planar = [circle]
        elif np.max(np.abs(point[:2])) >= a:
            planar = [point[:2].copy()]
        else:
            planar = []
            for axis in (0, 1):
                for sign in (-1., 1.):
                    projected = point[:2].copy()
                    projected[axis] = sign*a
                    planar.append(projected)
        candidates = []
        ambiguity = []
        for z in (0., H):
            options = np.array([np.r_[q, z] for q in planar])
            delta = point-options
            d2 = np.einsum('ij,ij->i', delta, delta)
            best = int(np.argmin(d2))
            # Compare squared distances using local cancellation bounds.
            tied = []
            for k, q in enumerate(options):
                terms = (q-options[best])*((q-point)+(options[best]-point))
                if abs(float(np.sum(terms))) <= 64*np.finfo(float).eps*np.sum(np.abs(terms)):
                    tied.append(k)
            distinct = len({tuple(options[k]) for k in tied}) > 1
            candidates.append(options[best]); ambiguity.append(distinct)
        for axis in (0, 1):
            for sign in (-1., 1.):
                q = point.copy()
                q[axis] = sign*a
                q[1-axis] = np.clip(point[1-axis], -a, a)
                q[2] = np.clip(point[2], 0., H)
                candidates.append(q); ambiguity.append(False)
        candidates.append(np.r_[circle, np.clip(point[2], 0., H)])
        ambiguity.append(radial == 0.)
        closest = np.asarray(candidates)
        delta = point-closest
        distances = np.linalg.norm(delta, axis=1)
        active = distances < self.particle_radius
        if np.any(distances == 0.):
            raise ContactGeometryError('particle centre on tool boundary; distance gradient undefined')
        if reject_active_ties and np.any(active & np.asarray(ambiguity)):
            names = [FEATURES[k] for k in np.flatnonzero(active & np.asarray(ambiguity))]
            raise ContactGeometryError('active multiple closest points: '+','.join(names))
        normals = delta/distances[:, None]
        return distances, normals, closest

    def geometry(self, positions, pose):
        """Endpoint geometry and conservative force, no friction/history update."""
        x = np.asarray(positions, float)
        if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
            raise ValueError('finite (N,3) positions required')
        origin, rotation = _pose(pose)
        local = (x-origin)@rotation
        gaps = np.full((len(x), 8), np.inf)
        normal_body = np.zeros((len(x), 8, 3))
        closest_body = np.zeros_like(normal_body)
        for i, point in enumerate(local):
            gaps[i, :7], normal_body[i, :7], closest_body[i, :7] = self._patches(point)
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
        normal_world[:, :7] = normal_body[:, :7]@rotation.T
        closest_world = closest_body.copy()
        closest_world[:, :7] = closest_body[:, :7]@rotation.T+origin
        return dict(local_positions_m=local, gaps_m=gaps, active=active,
                    normal_body=normal_body, normal_world=normal_world,
                    closest_points_m=closest_world, magnitude_N=magnitude,
                    force_by_feature_N=magnitude[:, :, None]*normal_world,
                    potential_J=.5*self.kn*overlap**2,
                    overlap_m=overlap)

    def _sweep(self, x0, x1, origin0, origin1, rotation0, rotation1, geometry0, geometry1):
        """Certify linear x/origin and geodesic tool rotation outside the solid.

        Distance to a fixed closed boundary is 1-Lipschitz. The body-coordinate
        curve speed is bounded by |dx-dorigin|+angle*max|x-origin|. A segment
        passes only when a sampled distance exceeds its whole interval bound;
        otherwise bisect. This cannot silently jump through a thin tool.
        """
        angle, rotation_at = _rotation_path(rotation0, rotation1)
        maximum_depth = 24
        maximum_queries = 4096
        evaluations = 0
        dx = x1-x0; do = origin1-origin0
        for i in range(len(x0)):
            speed_bound = (float(np.linalg.norm(dx[i]-do)) + angle*max(
                float(np.linalg.norm(x0[i]-origin0)), float(np.linalg.norm(x1[i]-origin1))))
            if speed_bound == 0.:
                continue
            start_distance = float(np.min(geometry0['gaps_m'][i, :7]))
            end_distance = float(np.min(geometry1['gaps_m'][i, :7]))
            stack = [(0., 1., start_distance, end_distance, 0)]
            while stack:
                lo, hi, dlo, dhi, depth = stack.pop()
                if max(dlo, dhi) > speed_bound*(hi-lo):
                    continue
                if depth >= maximum_depth:
                    raise ContactGeometryError('unresolved swept tool approach; subdivide process increment')
                middle = .5*(lo+hi)
                if evaluations >= maximum_queries:
                    raise ContactGeometryError('swept tool query budget exceeded; subdivide process increment')
                q = rotation_at(middle).T@(x0[i]+middle*dx[i]-origin0-middle*do)
                distances, _, _ = self._patches(q, reject_active_ties=False)
                evaluations += 1
                dm = float(np.min(distances))
                if dm > .5*speed_bound*(hi-lo):
                    continue
                stack.append((lo, middle, dlo, dm, depth+1))
                stack.append((middle, hi, dm, dhi, depth+1))
        return dict(interpolation='LINEAR_PARTICLE_AND_ORIGIN_SHORTEST_SO3_TOOL_ROTATION',
                    certified=True, intermediate_geometry_queries=evaluations,
                    maximum_subdivision_depth=maximum_depth,
                    maximum_geometry_queries=maximum_queries, tool_rotation_angle_rad=angle)

    def evaluate(self, x0, x1, pose0, pose1, dt, history=None, check_path=True):
        """Pure candidate contact force/history and unmodified signed work ledger.

        Finite wall displacement is defined at the old material body coordinate.
        Thus dx = wall_dx + R1*(y1-y0), exactly. Endpoint force quadrature is
        explicitly an approximation. ``normal_endpoint_work_error_J`` is its
        raw defect, never deposited as heat. Tangential projection/release
        energy is similarly an unassigned account, not physical friction heat.
        """
        x0 = np.asarray(x0, float); x1 = np.asarray(x1, float)
        if (x0.ndim != 2 or x0.shape[1] != 3 or x1.shape != x0.shape or
                not np.isfinite(x0).all() or not np.isfinite(x1).all() or
                not np.isfinite(dt) or dt <= 0.):
            raise ValueError('finite equal (N,3) positions and positive dt required')
        origin0, rotation0 = _pose(pose0); origin1, rotation1 = _pose(pose1)
        old = np.zeros((len(x0), 8, 3)) if history is None else np.asarray(history, float)
        if old.shape != (len(x0), 8, 3) or not np.isfinite(old).all():
            raise ValueError('finite (N,8,3) history required')
        g0 = self.geometry(x0, pose0); g1 = self.geometry(x1, pose1)
        sweep = (self._sweep(x0, x1, origin0, origin1, rotation0, rotation1, g0, g1)
                 if check_path else dict(certified=False, status='NOT_RUN_TANGENT_PERTURBATION_ONLY'))
        relative = np.empty_like(old)
        relative[:, :7] = (g1['local_positions_m']-g0['local_positions_m'])[:, None, :]
        relative[:, 7] = x1-x0
        normals = g1['normal_body']
        tangential_step = relative-np.sum(relative*normals, axis=2)[:, :, None]*normals
        projected_old = old-np.sum(old*normals, axis=2)[:, :, None]*normals
        projected_old = np.where(g1['active'][:, :, None], projected_old, 0.)
        trial = projected_old+tangential_step
        trial = np.where(g1['active'][:, :, None], trial, 0.)
        trial_norm = np.linalg.norm(trial, axis=2)
        limit = self.mu*g1['magnitude_N']/self.kt
        sliding = trial_norm > limit
        factor = np.ones_like(trial_norm)
        factor[sliding] = limit[sliding]/trial_norm[sliding]
        updated = trial*factor[:, :, None]
        plastic_step = trial-updated
        tangential_body = -self.kt*updated
        tangential_world = tangential_body.copy()
        tangential_world[:, :7] = tangential_body[:, :7]@rotation1.T
        friction = self.kt*np.sum(updated*plastic_step, axis=2)
        algorithmic = .5*self.kt*np.sum((updated-projected_old)**2, axis=2)
        energy0 = .5*self.kt*np.sum(old**2, axis=2)
        energy1 = .5*self.kt*np.sum(updated**2, axis=2)
        release = energy0-.5*self.kt*np.sum(projected_old**2, axis=2)
        tangential_relative_work = np.sum(tangential_body*relative, axis=2)
        tangential_balance = (tangential_relative_work+energy1-energy0+
                              friction+algorithmic+release)
        scale = (np.abs(tangential_relative_work)+energy0+energy1+
                 np.abs(friction)+algorithmic+np.abs(release))
        tolerance = 256*np.finfo(float).eps*scale
        if np.any(friction < -tolerance) or np.any(release < -tolerance):
            raise ArithmeticError('negative raw friction/release energy; candidate rejected')
        if np.any(np.abs(tangential_balance) > tolerance):
            raise ArithmeticError('tangential spring work identity failed')
        relative_world = relative.copy()
        relative_world[:, :7] = relative[:, :7]@rotation1.T
        wall_step = (x1-x0)[:, None, :]-relative_world
        normal = g1['force_by_feature_N']
        total = normal+tangential_world
        normal_relative_work = np.sum(normal*relative_world, axis=2)
        normal_error = g1['potential_J']-g0['potential_J']+normal_relative_work
        normal_particle_work = np.sum(normal*(x1-x0)[:, None, :], axis=2)
        normal_wall_work = np.sum(normal*wall_step, axis=2)
        particle_work = np.sum(total*(x1-x0)[:, None, :], axis=2)
        wall_work = np.sum(total*wall_step, axis=2)
        return dict(classification=CLASSIFICATION, features=FEATURES,
                    normal_model='SUM_OF_FINITE_PATCH_POTENTIALS_DUPLICATE_LIP_CONTRIBUTIONS_RETAINED',
                    friction_model='BODY_COORDINATE_TANGENTIAL_PENALTY_SPRING_COULOMB_RADIAL_RETURN',
                    work_quadrature='ENDPOINT_FORCE_FINITE_BODY_COORDINATE_WALL_DISPLACEMENT',
                    force_N=np.sum(total, axis=1), force_by_feature_N=total,
                    normal_force_by_feature_N=normal, tangential_force_by_feature_N=tangential_world,
                    history=updated, active=g1['active'], sliding=sliding & g1['active'],
                    gaps_m=g1['gaps_m'], normals=g1['normal_world'], normal_load_N=g1['magnitude_N'],
                    closest_points_m=g1['closest_points_m'], overlap_m=g1['overlap_m'],
                    penetration_m=g1['overlap_m'],
                    normal_potential0_J=g0['potential_J'], normal_potential1_J=g1['potential_J'],
                    tangential_potential0_J=energy0, tangential_potential1_J=energy1,
                    friction_dissipation_J=friction, algorithmic_dissipation_J=algorithmic,
                    release_energy_J=release, tangential_work_identity_error_J=tangential_balance,
                    normal_endpoint_work_error_J=normal_error,
                    normal_particle_work_J=normal_particle_work, normal_wall_work_J=normal_wall_work,
                    particle_work_J=particle_work, wall_work_J=wall_work,
                    relative_step_body_m=relative, wall_step_m=wall_step,
                    particle_reaction_on_tool_N=-np.sum(total[:, :7], axis=(0, 1)),
                    particle_reaction_torque_on_tool_Nm=-np.sum(np.cross(
                        (x1-origin1)[:, None, :], total[:, :7]), axis=(0, 1)),
                    particle_reaction_on_base_N=-np.sum(total[:, 7], axis=0),
                    sweep=sweep, dt_s=float(dt),
                    release_energy_policy='UNASSIGNED_ACCOUNT_NOT_AUTOMATICALLY_HEAT',
                    algorithmic_dissipation_policy='REPORTED_NOT_AUTOMATICALLY_HEAT')
