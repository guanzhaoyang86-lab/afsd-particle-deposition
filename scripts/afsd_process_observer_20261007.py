"""Passive rotation/torque/spreading observations for the new AFSD program.

All quantities use actual candidate/accepted data supplied by the caller.
No material velocity, graph, contact force, particle ID or history is changed.
Geometry-only material retention is never labelled metallurgical bonding.
"""
from pathlib import Path
import math
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:
    sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.interval_acceptance import validate_bore_crossing_state

CLASSIFICATION = 'PASSIVE_AFSD_ROTATION_TORQUE_AND_RETENTION_OBSERVER_20261007'


def observe_process(pa, state, result, motion, parameters):
    """Return JSON-friendly observations; no physical validity/acceptance claim.

    ``pa`` contains start-of-step particles; ``result`` is the actual evaluate
    result; ``state`` contains END-of-step bore-observation history. Tool poses
    and unwrapped angular displacement come from the common machine interval.
    Rotation is mass-weighted over component-1, material-labelled nonactuator
    particles, excluding only unresolved near-axis angular divisions.
    """
    x = np.asarray(result['position'], float)
    u = np.asarray(result['displacement'], float)
    ids, component = np.asarray(pa.particle_id), np.asarray(pa.component_id)
    mass = np.asarray(pa.m, float)
    actuator = np.asarray(pa.actuator).astype(bool)
    n = len(ids)
    if (x.shape != (n, 3) or u.shape != x.shape or mass.shape != (n,)
            or component.shape != (n,) or actuator.shape != (n,)
            or not np.isfinite(x).all() or not np.isfinite(u).all()
            or not np.isfinite(mass).all() or np.any(mass <= 0)
            or ids.dtype.kind not in 'iu' or len(np.unique(ids)) != n):
        raise ValueError('finite matching particles with persistent IDs and positive mass required')
    dt = float(motion.dt_s)
    t0, t1 = float(motion.t0), float(motion.t1)
    time_roundoff = 64*np.finfo(float).eps*max(abs(t0), abs(t1), abs(dt))
    if (not np.isfinite([dt, t0, t1]).all() or dt <= 0 or t1 <= t0
            or abs((t1-t0)-dt) > time_roundoff):
        raise ValueError('positive consistent finite motion interval required')
    origin0 = np.asarray(motion.pose(t0)['origin'], float)
    pose1 = motion.pose(t1)
    origin1, rotation1 = np.asarray(pose1['origin'], float), np.asarray(pose1['rotation'], float)
    if (origin0.shape != (3,) or origin1.shape != (3,) or rotation1.shape != (3, 3)
            or not np.isfinite(origin0).all() or not np.isfinite(origin1).all()
            or not np.isfinite(rotation1).all()
            or not np.allclose(rotation1.T@rotation1, np.eye(3), rtol=0., atol=1e-12)
            or abs(float(np.linalg.det(rotation1))-1.) > 1e-12
            or not np.allclose(rotation1[:, 2], [0., 0., 1.], rtol=0., atol=1e-12)):
        raise ValueError('finite vertical-axis tool poses required')
    radius, bore, dx = (float(parameters[k]) for k in ('tool_radius', 'bore_width', 'dx'))
    substrate_top = float(parameters.get('initial_substrate_top_z_m', 0.))
    if not np.isfinite([radius, bore, dx, substrate_top]).all() or min(radius, bore, dx) <= 0:
        raise ValueError('finite positive particle/tool dimensions required')
    description = motion.to_dict()
    if 'average_omega_rad_s' in description:
        tool_omega = float(description['average_omega_rad_s'])
    else:
        tool_omega = (float(description['unwrapped_angle1_rad'])-
                      float(description['unwrapped_angle0_rad']))/dt
    if not math.isfinite(tool_omega):
        raise ValueError('finite unwrapped tool angular velocity required')
    history = state['bore_observation']
    validate_bore_crossing_state(history)
    order = np.argsort(ids)
    if (ids[order].astype(int).tolist() != history['particle_ids']
            or component[order].astype(int).tolist() != history['component_ids']
            or mass[order].tolist() != history['masses_kg']
            or not math.isclose(float(history['last_time_s']), t1, rel_tol=0., abs_tol=1e-12)
            or history['bore_width_m'] != bore):
        raise ValueError('end-step observation history must match persistent particles and motion')

    eligible = np.isin(ids, history['eligible_ids']) & (component == 1)
    crossed = eligible & np.isin(ids, history['crossing_ids'])
    nonactuator_bar = (component == 1) & ~actuator
    r = x-origin1
    radial_squared = np.sum(r[:, :2]**2, axis=1)
    # Roundoff-scale axis exclusion is reported, not used to move particles.
    axis_tolerance = 64*np.finfo(float).eps*max(radius, bore, dx,
        float(np.max(np.abs(x), initial=0.)), float(np.max(np.abs(origin1), initial=0.)))
    rotating = nonactuator_bar & (radial_squared > axis_tolerance**2)
    relative_velocity = u/dt-(origin1-origin0)/dt
    cross_z = r[:, 0]*relative_velocity[:, 1]-r[:, 1]*relative_velocity[:, 0]
    omega_values = cross_z[rotating]/radial_squared[rotating]
    omega_mass = math.fsum(mass[rotating].tolist())
    material_omega = (float(np.sum(mass[rotating]*omega_values))/omega_mass
                      if omega_mass else None)
    inertia = float(np.sum(mass[rotating]*radial_squared[rotating]))
    inertia_omega = float(np.sum(mass[rotating]*cross_z[rotating]))/inertia if inertia else None

    local = r@rotation1
    outside_square = np.max(np.abs(local[:, :2]), axis=1) >= bore/2
    above_substrate_below_shoulder = (x[:, 2] > substrate_top) & (local[:, 2] < 0.)
    under_annular_shoulder = crossed & outside_square & (radial_squared <= radius**2) & above_substrate_below_shoulder
    entirely_clear_of_shoulder_projection = (radial_squared > (radius+dx/2)**2)
    retained_candidate = crossed & entirely_clear_of_shoulder_projection & above_substrate_below_shoulder

    def population(mask):
        return dict(particle_count=int(mask.sum()), particle_ids=ids[mask].astype(int).tolist(),
                    mass_kg=math.fsum(mass[mask].tolist()))

    contact = result.get('contact', {})
    force = contact.get('force_by_feature_N')
    if force is None:
        torque = dict(status='NOT_AVAILABLE_NO_SAVED_FEATURE_FORCES',
                      force_reconstructed=False)
    else:
        force = np.asarray(force, float)
        if force.shape != (n, 8, 3) or not np.isfinite(force).all():
            raise ValueError('finite saved contact feature forces with shape N,8,3 required')
        bore_force = force[:, 2:6, :].sum(axis=1)
        particle_torque = r[:, 0]*bore_force[:, 1]-r[:, 1]*bore_force[:, 0]
        all_torque = math.fsum(particle_torque.tolist())
        torque = dict(status='OBSERVED_SAVED_FEATURE_FORCES', feature_indices=[2, 3, 4, 5],
            on_all_deformable_material_z_Nm=all_torque,
            reaction_on_tool_z_Nm=-all_torque,
            on_nonactuator_bar_z_Nm=math.fsum(particle_torque[nonactuator_bar].tolist()),
            on_actuator_bar_z_Nm=math.fsum(particle_torque[(component == 1) & actuator].tolist()),
            lever_arm='ENDPOINT_PARTICLE_POSITION_MINUS_ENDPOINT_TOOL_ORIGIN',
            discrete_contact_work_identity_certified=False, force_reconstructed=False)
    report = dict(classification=CLASSIFICATION, status='PASSIVE_OBSERVATION_ONLY',
        time0_s=t0, time1_s=t1, physical_solver_evaluations=0, material_trajectory_prescribed=False,
        acceptance_changed=False,
        rotation=dict(tool_average_omega_rad_s=tool_omega,
            mass_weighted_material_omega_rad_s=material_omega,
            inertia_weighted_material_omega_rad_s=inertia_omega,
            mass_weighted_material_to_tool_ratio=(material_omega/tool_omega
                if material_omega is not None and tool_omega != 0. else None),
            tool_minus_mass_weighted_material_rad_s=(tool_omega-material_omega
                if material_omega is not None else None),
            population=population(nonactuator_bar), resolved_population=population(rotating),
            axis_excluded_population=population(nonactuator_bar & ~rotating),
            axis_division_tolerance_m=axis_tolerance,
            tool_translation_velocity_m_s=((origin1-origin0)/dt).tolist(),
            estimator='ENDPOINT_RADIUS_AND_TRANSLATION_SUBTRACTED_STEP_AVERAGE_CHORD_VELOCITY',
            imposed_material_rotation=False),
        bore_wall_torque=torque,
        supplied_material=dict(eligible=population(eligible), eligible_crossed=population(crossed),
            eligible_not_yet_crossed=population(eligible & ~crossed),
            crossed_inside_square_projection=population(crossed & ~outside_square),
            crossed_under_annular_shoulder=population(under_annular_shoulder),
            geometric_retention_candidate=population(retained_candidate),
            geometric_retention_is_metallurgical_bonding=False,
            initial_substrate_top_z_m=substrate_top,
            initial_substrate_top_source=('EXPLICIT_PARAMETER' if 'initial_substrate_top_z_m' in parameters
                                          else 'EXISTING_GEOMETRY_DEFAULT_Z_ZERO'),
            complete_shoulder_clearance_radius_m=radius+dx/2,
            masks_use='ENDPOINT_CENTRES_EXCEPT_RADIUS_PLUS_HALF_DX_FOR_FULL_HORIZONTAL_CLEARANCE'),
        limitations=[
            'Finite-step chord angular velocity is an observer estimate, not a prescribed spin or continuum convergence claim.',
            'Near-axis particles have undefined resolved angular division and are explicitly excluded from angular averages.',
            'Only initially eligible, historically crossed persistent IDs count as supplied material; initial lower bar is excluded.',
            'Geometric retention candidates do not establish stable adhesion, metallurgical bonding or experimental deposition.',
            'Endpoint feature-force torque is diagnostic and is not substituted for the discrete contact energy accounting.',
        ])
    return report
