"""Explicit centre-collocated rigid-contact companion; no force-law changes.

Face projection is diagnostic geometry, NOT a complete closest-solid-point search.
No particle spin DOF and no force transfer from a reconstructed face velocity.
SI only. Tracer/AS/RS labels are deliberately absent from all interfaces.
"""
import numpy as np
from .geometry import surfaces, assert_outside_tool, tool_pose

FEATURE_NAMES = ('BASE', 'SHOULDER', 'BORE_X_MINUS', 'BORE_X_PLUS',
                 'BORE_Y_MINUS', 'BORE_Y_PLUS')
CLASSIFICATION = 'CENTRE_COLLOCATED_APPROXIMATE_FACE_POINT_DIAGNOSTIC'


def array(value, shape, name):
    a = np.asarray(value, dtype=float)
    if a.shape != shape or not np.isfinite(a).all():
        raise ValueError('invalid/nonfinite '+name)
    return a


def rigid_velocity(points, origin, translation, omega):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('invalid/nonfinite points')
    origin = array(origin, (3,), 'origin')
    translation = array(translation, (3,), 'translation')
    omega = array(omega, (3,), 'omega')
    return translation+np.cross(omega, points-origin)


def contact_rows(x, particle_id, t, p):
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
        raise ValueError('invalid/nonfinite particle position')
    ids = np.asarray(particle_id)
    if (ids.shape != (len(x),) or ids.dtype.kind not in 'iu' or
            len(np.unique(ids)) != len(ids)):
        raise ValueError('unique integer persistent particle IDs required')
    keys = ('dx', 'traverse', 'omega', 'layer', 'tool_top', 'tool_radius',
            'bore_width', 'substrate_depth', 'base_temperature', 'tool_temperature')
    if not np.isfinite([t]+[p[k] for k in keys]).all():
        raise ValueError('nonfinite time/geometry parameter')
    if (p['dx'] <= 0 or p['bore_width'] <= 0 or p['substrate_depth'] <= 0 or
            p['tool_top'] <= p['layer'] or p['tool_radius'] <= p['bore_width']/np.sqrt(2)):
        raise ValueError('invalid finite open-square tool geometry')
    assert_outside_tool(x, t, p)
    origin, _ = tool_pose(t, p['traverse'], p['omega'], p['layer'])
    index, feature, gap, normal, application, face, wall, facewall = ([] for _ in range(8))
    origins, translations, omegas = [], [], []
    for fid, (indices, gaps, normals, old_wall_velocity, _) in enumerate(surfaces(x, t, p)):
        o = np.array([0., 0., -p['substrate_depth']]) if fid == 0 else origin
        tr = np.zeros(3) if fid == 0 else np.array([p['traverse'], 0., 0.])
        om = np.zeros(3) if fid == 0 else np.array([0., 0., p['omega']])
        centre = x[indices]
        projected = centre-gaps[:, None]*normals
        vw = rigid_velocity(centre, o, tr, om)
        vc = rigid_velocity(projected, o, tr, om)
        if not np.allclose(vw, old_wall_velocity, rtol=1e-13, atol=1e-14):
            raise ValueError('old centre-collocated kinematics changed')
        for j, particle in enumerate(indices):
            index.append(int(particle)); feature.append(fid); gap.append(gaps[j])
            normal.append(normals[j]); application.append(centre[j]); face.append(projected[j])
            wall.append(vw[j]); facewall.append(vc[j])
            origins.append(o); translations.append(tr); omegas.append(om)
    index = np.asarray(index, dtype=np.int64)
    feature = np.asarray(feature, dtype=np.int64)
    gaps = np.asarray(gap)
    # This is geometry overlap, NOT clipping a force, work, temperature or state.
    penetration = p['dx']/2-gaps
    if np.any(penetration < 0):
        raise ValueError('broad-phase returned a separated contact')
    triples = lambda a: np.asarray(a, dtype=float).reshape(-1, 3)
    return dict(classification=CLASSIFICATION, particle_index=index,
        particle_id=ids[index].copy(), feature_id=feature,
        contact_id=np.array([f'{int(ids[i])}:{int(f)}' for i, f in zip(index, feature)], dtype=str),
        normal=triples(normal), gap=gaps, penetration=penetration,
        application_point=triples(application), contact_point=triples(face),
        wall_velocity=triples(wall), wall_velocity_at_face=triples(facewall),
        body_origin=triples(origins), translation_velocity=triples(translations), omega=triples(omegas))


def relative_kinematics(rows, velocity, velocity_gradient=None):
    velocity = np.asarray(velocity, dtype=float)
    if velocity.ndim != 2 or velocity.shape[1] != 3 or not np.isfinite(velocity).all():
        raise ValueError('invalid/nonfinite material velocity')
    index = rows['particle_index']
    if len(index) and (np.min(index) < 0 or np.max(index) >= len(velocity)):
        raise ValueError('contact index outside particle inventory')
    material = velocity[index]
    n = rows['normal']
    relative = material-rows['wall_velocity']
    tangent = relative-np.sum(relative*n, axis=1)[:, None]*n
    material_face = relative_face = None
    if velocity_gradient is not None:
        gradient = array(velocity_gradient, (len(velocity), 3, 3), 'material velocity gradient')
        offset = rows['contact_point']-rows['application_point']
        material_face = material+np.einsum('nij,nj->ni', gradient[index], offset)
        relative_face = material_face-rows['wall_velocity_at_face']
    return dict(material_velocity=material, relative_velocity=relative,
                tangential_relative_velocity=tangent, material_velocity_at_face=material_face,
                relative_velocity_at_face=relative_face,
                face_reconstruction_status='DIAGNOSTIC_ONLY_NOT_FORCE_TRANSFER' if velocity_gradient is not None
                                           else 'NOT_SOLVED_NO_VELOCITY_GRADIENT')


def contact_power(rows, forces, velocity):
    f = array(forces, rows['normal'].shape, 'contact forces')
    k = relative_kinematics(rows, velocity)
    n = rows['normal']
    fn = np.sum(f*n, axis=1)[:, None]*n
    ft = f-fn
    raw = -np.sum(ft*k['tangential_relative_velocity'], axis=1)
    # Normal relative work may store/release penalty energy; never call it heat.
    normal_relative = -np.sum(fn*k['relative_velocity'], axis=1)
    ww = np.sum(f*rows['wall_velocity'], axis=1)
    wm = np.sum(f*k['material_velocity'], axis=1)
    delta_torque = np.cross(rows['contact_point']-rows['application_point'], f)
    delta_power = np.sum(f*(rows['wall_velocity_at_face']-rows['wall_velocity']), axis=1)
    torque_power = np.sum(rows['omega']*delta_torque, axis=1)
    return dict(friction_power_raw_W=raw, normal_relative_power_W=normal_relative,
                wall_power_W=ww, material_power_W=wm,
                closure_residual_W=ww-wm-raw-normal_relative,
                face_minus_application_torque_Nm=delta_torque,
                face_minus_application_power_W=delta_power,
                torque_power_difference_W=delta_power-torque_power)
