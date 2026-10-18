"""Exact finite open-square boundary geometry; independent approximate penalty.

No solver integration. No friction, thermal or particle motion. Sharp features
are not rounded. Active ambiguous distance gradients are explicitly rejected.
The one-body distance potential is NOT the old sum of per-face springs.
"""
import numpy as np


CLASSIFICATION = 'FINITE_OPEN_SQUARE_DISTANCE_APPROXIMATE_COMPONENT_NOT_AFSD'


class NonsmoothContactError(ValueError):
    def __init__(self, query):
        super().__init__('active NONSMOOTH finite-tool distance: '+query['normal_status'])
        self.query = query


def _vector(value, name):
    result = np.asarray(value, dtype=float)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError('finite three-vector required: '+name)
    return result


class FiniteOpenSquareTool:
    """Local solid: finite circular cylinder minus a through square bore.

    Origin is centre of the lower shoulder. Rotation maps body to laboratory
    coordinates; geometry does not prescribe any material motion.
    """
    def __init__(self, radius, bore_width, height):
        if (not np.isfinite([radius, bore_width, height]).all() or
                radius <= 0. or bore_width <= 0. or height <= 0. or
                radius <= bore_width/np.sqrt(2.)):
            raise ValueError('positive finite dimensions and circle outside square required')
        self.radius = float(radius)
        self.half_bore = .5*float(bore_width)
        self.height = float(height)

    def _candidates(self, point):
        """Exact closest points to each finite boundary component, with ties.

        Clamps here are geometric projections onto closed segments, not state
        or dissipation clipping. The circle contains the entire square bore.
        """
        a, R, H = self.half_bore, self.radius, self.height
        xy = point[:2]
        radial = np.linalg.norm(xy)
        circle = xy*(R/radial) if radial > 0. else np.array([R, 0.])
        if radial > R:
            planar = [circle]
        elif np.max(np.abs(xy)) >= a:
            planar = [xy.copy()]
        else:
            planar = []
            for axis in (0, 1):
                for sign in (-1, 1):
                    projected = xy.copy()
                    projected[axis] = sign*a
                    planar.append(projected)
        result = []
        for name, z, nz in [('SHOULDER', 0., -1.), ('TOP', H, 1.)]:
            for projected in planar:
                result.append((name, np.r_[projected, z], np.array([0., 0., nz])))
        for axis, letter in ((0, 'X'), (1, 'Y')):
            for sign, side in ((-1, 'MINUS'), (1, 'PLUS')):
                projected = point.copy()
                projected[axis] = sign*a
                projected[1-axis] = np.clip(point[1-axis], -a, a)
                projected[2] = np.clip(point[2], 0., H)
                normal = np.zeros(3); normal[axis] = -sign
                result.append(('BORE_'+letter+'_'+side, projected, normal))
        result.append(('OUTER_CYLINDER', np.r_[circle, np.clip(point[2], 0., H)],
                       np.r_[circle/R, 0.]))
        return result

    def query(self, point, origin=(0., 0., 0.), rotation=None,
              translation_velocity=(0., 0., 0.), angular_velocity=(0., 0., 0.)):
        point = _vector(point, 'point')
        origin = _vector(origin, 'origin')
        translation = _vector(translation_velocity, 'translation velocity')
        omega = _vector(angular_velocity, 'angular velocity')
        rotation = np.eye(3) if rotation is None else np.asarray(rotation, dtype=float)
        if (rotation.shape != (3, 3) or not np.isfinite(rotation).all() or
                not np.allclose(rotation.T@rotation, np.eye(3), atol=1e-12, rtol=0.) or
                abs(np.linalg.det(rotation)-1.) > 1e-12):
            raise ValueError('proper orthonormal three-dimensional rotation required')
        local = rotation.T@(point-origin)
        candidates = self._candidates(local)
        distances = np.array([np.linalg.norm(local-q) for _, q, _ in candidates])
        distance = float(np.min(distances))
        tol = 64*np.finfo(float).eps*max(self.radius, self.height, np.linalg.norm(local))
        # Compare squared distances as (qa-qb).((qa-x)+(qb-x)). Comparing
        # nearly equal norms with a global absolute tolerance falsely makes a
        # smooth point 1e-8 m from a lip look like a two-point medial axis.
        # The dot-product rounding bound is local to the comparison terms.
        selected = []
        for entry in candidates:
            qa = entry[1]
            dominated = False
            for _, qb, _ in candidates:
                terms = (qa-qb)*((qa-local)+(qb-local))
                roundoff = 64*np.finfo(float).eps*np.sum(np.abs(terms))
                if np.sum(terms) > roundoff:
                    dominated = True
                    break
            if not dominated:
                selected.append(entry)
        if not selected:
            raise ArithmeticError('inconsistent finite-tool closest-point comparisons')
        distinct = []
        for _, closest, _ in selected:
            if not any(np.linalg.norm(closest-other) <= tol for other in distinct):
                distinct.append(closest)
        inside = (0. < local[2] < self.height and
                  np.linalg.norm(local[:2]) < self.radius and
                  np.max(np.abs(local[:2])) > self.half_bore)
        signed = -distance if inside else distance
        normal = None
        if len(distinct) > 1:
            status = 'NONSMOOTH_MULTIPLE_CLOSEST_POINTS'
        elif distance > tol:
            direction = (local-distinct[0])/distance
            normal = rotation@(-direction if inside else direction)
            status = 'UNIQUE_DISTANCE_GRADIENT'
        else:
            # On a sharp boundary the adjacent face normals disagree. Do not
            # pick the first face or invent a smoothed average normal.
            normals = [entry[2] for entry in selected]
            if all(np.linalg.norm(n-normals[0]) <= 1e-12 for n in normals):
                normal = rotation@normals[0]
                status = 'UNIQUE_SMOOTH_BOUNDARY_NORMAL'
            else:
                status = 'NONSMOOTH_BOUNDARY_EDGE_OR_CORNER'
        closest = np.array([origin+rotation@q for q in distinct])
        velocity_at_closest = translation+np.cross(omega, closest-origin)
        velocity_at_application = translation+np.cross(omega, point-origin)
        return dict(classification=CLASSIFICATION, body_id='FINITE_TOOL',
                    signed_distance_m=signed, inside_solid=bool(inside),
                    boundary_features=tuple(sorted(set(entry[0] for entry in selected))),
                    closest_points_m=closest, normal=normal, normal_status=status,
                    application_point_m=point.copy(), body_origin_m=origin.copy(),
                    translation_velocity_m_s=translation.copy(), angular_velocity_rad_s=omega.copy(),
                    wall_velocity_at_closest_m_s=velocity_at_closest,
                    wall_velocity_at_application_m_s=velocity_at_application,
                    boundary_is_sharp=len(set(entry[0] for entry in selected)) > 1,
                    closest_point_merge_tolerance_m=tol,
                    squared_distance_comparison_roundoff_factor=64)


def normal_penalty(query, particle_radius, stiffness):
    """Conservative one-body distance penalty. Units: m, N/m, N and J.

    This is a separate named variant, not a backwards-compatible contact_rows
    replacement. A component query is not evidence of time-integrator work.
    """
    if (query.get('classification') != CLASSIFICATION or
            not np.isfinite([particle_radius, stiffness]).all() or
            particle_radius <= 0. or stiffness <= 0.):
        raise ValueError('component query and positive finite radius/stiffness required')
    signed = float(query['signed_distance_m'])
    if not np.isfinite(signed):
        raise ValueError('nonfinite signed distance')
    overlap = particle_radius-signed
    active = overlap > 0.
    if active and query['normal'] is None:
        raise NonsmoothContactError(query)
    delta = overlap if active else 0.
    force = stiffness*delta*query['normal'] if active else np.zeros(3)
    return dict(classification=CLASSIFICATION, active=active,
                signed_penetration_m=overlap, penetration_m=delta,
                normal_force_N=force, potential_J=.5*stiffness*delta*delta,
                model='SINGLE_BODY_DISTANCE_NOT_SUM_OF_FACE_SPRINGS')
