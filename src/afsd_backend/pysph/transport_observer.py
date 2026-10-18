"""Read-only, ID-aligned observations between two accepted material states.

Particle-centre crossings are not deposited volume, a 2-D area fraction, or a
continuum flux reconstruction. This module never advances a particle or creates
mass. Callers must establish acceptance of both input states before calling it.
"""
import hashlib
import json

import numpy as np

from .axial_actuator import control


CLASSIFICATION = 'ACCEPTED_ENDPOINT_TRANSPORT_OBSERVATION_NOT_DEPOSITION'


def _get(pa, key):
    value = np.asarray(getattr(pa, key))
    if value.ndim != 1 or not np.isfinite(value).all():
        raise ValueError('finite one-dimensional property required: '+key)
    return value


def _snapshot(pa):
    ids = _get(pa, 'particle_id')
    if len(ids) == 0 or len(np.unique(ids)) != len(ids) or not np.equal(ids, np.floor(ids)).all():
        raise ValueError('nonempty unique integer particle IDs required')
    order = np.argsort(ids)
    keys = ('particle_id', 'component_id', 'material_id', 'x0', 'y0', 'z0',
            'x', 'y', 'z', 'm', 'rho', 'ep')
    values = {key: _get(pa, key)[order].copy() for key in keys}
    if any(len(v) != len(ids) for v in values.values()):
        raise ValueError('particle property length mismatch')
    if np.any(values['m'] <= 0) or np.any(values['rho'] <= 0) or np.any(values['ep'] < 0):
        raise ValueError('positive mass/density and nonnegative plastic strain required')
    if not np.isin(values['component_id'], [0, 1]).all():
        raise ValueError('only deformable substrate=0 and bar=1 are registered')
    labels = {}
    for key in ('tracer_type', 'wholebar_tag', 'as_rs_tag', 'actuator'):
        if hasattr(pa, key):
            value = _get(pa, key)
            if len(value) != len(ids): raise ValueError('label length mismatch')
            labels[key] = value[order].tolist()
    values['labels'] = labels
    return values


def _selection(mask, ids, mass):
    return dict(particle_ids=ids[mask].tolist(), particle_count=int(mask.sum()),
                mass_kg=float(mass[mask].sum()))


def observe(before_pa, after_pa, t0, t1, p, cumulative=None):
    """Return detached JSON data; ``cumulative=result['cumulative']`` chains it.

    Persistent mass/identity and original coordinates must be unchanged. This
    finite-material observer deliberately rejects insertion/deletion; an inlet
    needs a separately specified source and crossing operator, not guessed mass.
    ``z >= layer`` is the upper half-space. A downward event ends strictly below
    it, and an upward event ends on or above it. Within-interval recrossings are
    NOT_RESOLVED by endpoint sampling.
    """
    t0, t1 = float(t0), float(t1)
    if not np.isfinite([t0, t1]).all() or t0 < 0 or t1 <= t0:
        raise ValueError('ordered finite nonnegative observation times required')
    dt = t1-t0
    layer, traverse, dx = (float(p[k]) for k in ('layer', 'traverse', 'dx'))
    if not np.isfinite([layer, traverse, dx]).all() or dx <= 0:
        raise ValueError('finite geometry and positive spacing required')
    a, b = _snapshot(before_pa), _snapshot(after_pa)
    for key in ('particle_id', 'component_id', 'material_id', 'x0', 'y0', 'z0', 'm'):
        if not np.array_equal(a[key], b[key]):
            raise ValueError('fixed material identity/mass changed: '+key)
    if np.any(b['ep'] < a['ep']):
        raise ValueError('accepted equivalent plastic strain decreased')
    ids, mass = a['particle_id'].astype(np.int64), a['m']
    xyz0 = np.column_stack([a[k] for k in 'xyz'])
    xyz1 = np.column_stack([b[k] for k in 'xyz'])
    lab = xyz1-xyz0
    origin0, origin1 = np.array([traverse*t0, 0., layer]), np.array([traverse*t1, 0., layer])
    relative = (xyz1-origin1)-(xyz0-origin0)
    bar = a['component_id'] == 1
    command = [control(a['z0'], t, 0., p['bar_length'], p['actuator_depth'],
                       p['feed'], p['release_time']) for t in (t0, t1)]
    active0, active1 = (c['active'] for c in command)
    free = ~(active0 | active1)
    down = bar & (a['z'] >= layer) & (b['z'] < layer)
    up = bar & (a['z'] < layer) & (b['z'] >= layer)
    delta_ep = b['ep']-a['ep']

    groups = {}
    for name, mask in (('bar', bar), ('substrate', ~bar), ('free_bar', bar & free),
                       ('free_substrate', ~bar & free), ('actuator_at_either_endpoint', ~free)):
        row = _selection(mask, ids, mass)
        total = float(mass[mask].sum())
        row.update(max_lab_displacement_m=float(np.linalg.norm(lab[mask], axis=1).max(initial=0.)),
                   max_tool_translation_displacement_m=float(np.linalg.norm(relative[mask], axis=1).max(initial=0.)),
                   lab_mass_weighted_displacement_m=(np.sum(mass[mask, None]*lab[mask], axis=0)/total).tolist() if total else None,
                   tool_mass_weighted_displacement_m=(np.sum(mass[mask, None]*relative[mask], axis=0)/total).tolist() if total else None,
                   volume_before_m3=float(np.sum(mass[mask]/a['rho'][mask])),
                   volume_after_m3=float(np.sum(mass[mask]/b['rho'][mask])))
        row['max_lab_displacement_over_dx'] = row['max_lab_displacement_m']/dx
        groups[name] = row

    signature = hashlib.sha256(json.dumps(dict(ids=ids.tolist(), mass=mass.tolist(),
        layer=layer, traverse=traverse), sort_keys=True, allow_nan=False).encode()).hexdigest()
    previous = dict(cumulative or {})
    if cumulative is not None:
        if previous.get('identity_plane_signature') != signature or previous.get('last_time_s') != t0:
            raise ValueError('cumulative history identity/plane or endpoint time mismatch')
    old_down, old_up = (set(previous.get(k, [])) for k in ('downward_unique_ids', 'upward_unique_ids'))
    if not (old_down | old_up).issubset(set(ids.tolist())):
        raise ValueError('unknown ID in cumulative history')
    new_down, new_up = old_down | set(ids[down].tolist()), old_up | set(ids[up].tolist())
    history = dict(identity_plane_signature=signature, first_time_s=previous.get('first_time_s', t0),
        last_time_s=t1, downward_unique_ids=sorted(new_down), upward_unique_ids=sorted(new_up),
        downward_unique_mass_kg=float(mass[np.isin(ids, list(new_down))].sum()),
        upward_unique_mass_kg=float(mass[np.isin(ids, list(new_up))].sum()),
        downward_event_count=int(previous.get('downward_event_count', 0)+down.sum()),
        upward_event_count=int(previous.get('upward_event_count', 0)+up.sum()),
        downward_event_mass_kg=float(previous.get('downward_event_mass_kg', 0.)+mass[down].sum()),
        upward_event_mass_kg=float(previous.get('upward_event_mass_kg', 0.)+mass[up].sum()),
        interpretation='ENDPOINT_CENTRE_EVENTS; REPEATED_CROSSINGS_COUNT_AGAIN_IN_EVENT_MASS; NO_EXTERNAL_SOURCE')
    crossing = dict(plane_z_m=layer, normal_lab=[0., 0., 1.],
        downward=_selection(down, ids, mass), upward=_selection(up, ids, mass),
        signed_downward_event_mass_kg=float(mass[down].sum()-mass[up].sum()),
        signed_downward_event_mass_rate_kg_s=float((mass[down].sum()-mass[up].sum())/dt),
        status='CENTER_CROSSING_OBSERVED_NOT_DEPOSITION' if np.any(down | up) else 'NO_CENTER_CROSSING_NOT_DEPOSITION',
        within_interval_recrossings='NOT_RESOLVED_ENDPOINT_SAMPLING',
        interpretation='PARTICLE_CENTRE_EVENT_MASS; NOT_CONTINUUM_FLUX_OR_DEPOSITED_VOLUME')
    result = dict(classification=CLASSIFICATION, input_acceptance='CALLER_MUST_VERIFY_BOTH_ACCEPTED',
        time_before_s=t0, time_after_s=t1, interval_s=dt, groups=groups, crossing=crossing,
        plasticity=dict(existing_before=_selection(a['ep'] > 0, ids, mass),
            nonzero_after=_selection(b['ep'] > 0, ids, mass),
            increased_this_interval=_selection(delta_ep > 0, ids, mass),
            continued_increase=_selection((a['ep'] > 0) & (delta_ep > 0), ids, mass)),
        per_particle=dict(particle_id=ids.tolist(), component_id=a['component_id'].tolist(),
            material_id=a['material_id'].tolist(), position_before_m=xyz0.tolist(), position_after_m=xyz1.tolist(),
            lab_displacement_m=lab.tolist(), tool_translation_displacement_m=relative.tolist(),
            mass_kg=mass.tolist(), volume_before_m3=(mass/a['rho']).tolist(), volume_after_m3=(mass/b['rho']).tolist(),
            equivalent_plastic_strain_before=a['ep'].tolist(), equivalent_plastic_strain_after=b['ep'].tolist(),
            equivalent_plastic_strain_increment=delta_ep.tolist(), free_at_both_endpoints=free.tolist(),
            actuator_active_before=active0.tolist(), actuator_active_after=active1.tolist(),
            observation_labels_before=a['labels'], observation_labels_after=b['labels']),
        tool_reference=dict(translation_only=True, rotates=False, origin_before_m=origin0.tolist(), origin_after_m=origin1.tolist()),
        cumulative=history, external_mass_source_kg=0.,
        macro_transport_status='NOT_ESTABLISHED_BY_ENDPOINT_MOTION',
        deposition_status='NOT_MEASURED', area_fraction_status='NOT_MEASURED_NO_3D_COUNT_TO_2D_AREA')
    # Guard the stated output interface, without rounding or clipping diagnostics.
    json.dumps(result, allow_nan=False)
    return result
