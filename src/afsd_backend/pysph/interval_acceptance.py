"""Pure interval energy gates and passive, geometric square-bore observation.

These checks neither solve mechanics nor certify a physical process.  The
energy contract is frozen here before new trials: each step and every prefix
must pass, with only ONE absolute allowance for the whole independent history.
The sum of absolute step residuals is also gated, so cancellation cannot make
an otherwise failed history acceptable.  No algorithmic loss is subtracted.
"""
from copy import deepcopy
import math
import numpy as np

ENERGY_CONTRACT = 'ONE_ABSOLUTE_ALLOWANCE_AND_NO_STEP_CANCELLATION_V1'
ENERGY_ABSOLUTE_J = 1e-7
ENERGY_RELATIVE = .05
CROSSING_CONTRACT = 'INITIAL_ABOVE_BAR_CENTRE_THROUGH_MOVING_OPEN_SQUARE_V1'


def _finite(value, name):
    try:
        value = float(value)
    except (TypeError, OverflowError) as exc:
        raise ValueError('invalid '+name) from exc
    if not math.isfinite(value):
        raise ValueError('nonfinite '+name)
    return value


def _energy_record(energy):
    if not {'global_raw_residual_J', 'absolute_work_scale_J'}.issubset(energy):
        raise ValueError('missing raw energy residual or absolute work scale')
    raw = _finite(energy['global_raw_residual_J'], 'raw energy residual')
    scale = _finite(energy['absolute_work_scale_J'], 'absolute work scale')
    if scale < 0:
        raise ValueError('negative absolute work scale')
    return dict(global_raw_residual_J=raw, absolute_work_scale_J=scale)


def _energy_summary(records):
    raw = math.fsum(r['global_raw_residual_J'] for r in records)
    absolute = math.fsum(abs(r['global_raw_residual_J']) for r in records)
    scale = math.fsum(r['absolute_work_scale_J'] for r in records)
    limit = ENERGY_ABSOLUTE_J+ENERGY_RELATIVE*scale
    if not all(math.isfinite(x) for x in (raw, absolute, scale, limit)):
        raise ValueError('nonfinite cumulative energy')
    return dict(steps=len(records), cumulative_raw_residual_J=raw,
                cumulative_absolute_step_residual_J=absolute,
                cumulative_absolute_work_scale_J=scale,
                interval_limit_J=limit)


def initial_interval_energy_state():
    """Call only at a separately declared independent initialization."""
    return dict(contract=ENERGY_CONTRACT, absolute_allowance_J=ENERGY_ABSOLUTE_J,
                relative_allowance=ENERGY_RELATIVE, records=[],
                **_energy_summary([]))


def _energy_gates(records, summary):
    last = records[-1]
    return dict(
        step=abs(last['global_raw_residual_J']) <=
             ENERGY_ABSOLUTE_J+ENERGY_RELATIVE*last['absolute_work_scale_J'],
        interval_signed=abs(summary['cumulative_raw_residual_J']) <=
                        summary['interval_limit_J'],
        interval_no_cancellation=summary['cumulative_absolute_step_residual_J'] <=
                                 summary['interval_limit_J'])


def validate_interval_energy_state(state):
    """Validate a serialized accepted history, including every accepted prefix.

Requiring the full two-scalar record history prevents a checkpoint from
presenting only a reset or inconsistent accumulator.  The caller must also
match ``steps`` to its physical step count and independent initialization.
"""
    required = set(initial_interval_energy_state())
    if set(state) != required:
        raise ValueError('interval energy state fields differ from contract')
    if (state['contract'] != ENERGY_CONTRACT or
            state['absolute_allowance_J'] != ENERGY_ABSOLUTE_J or
            state['relative_allowance'] != ENERGY_RELATIVE):
        raise ValueError('interval energy contract changed')
    if not isinstance(state['records'], list):
        raise ValueError('interval energy records must be a list')
    records = []
    for row in state['records']:
        if set(row) != {'global_raw_residual_J', 'absolute_work_scale_J'}:
            raise ValueError('invalid interval energy record fields')
        records.append(_energy_record(row))
        if not all(_energy_gates(records, _energy_summary(records)).values()):
            raise ValueError('checkpoint contains a failed interval energy prefix')
    summary = _energy_summary(records)
    for key, value in summary.items():
        if state[key] != value:
            raise ValueError('inconsistent interval energy field: '+key)
    return True


def assess_interval_energy(previous, energy):
    """Evaluate without changing ``previous``; rejection supplies no next state.

Malformed input raises ValueError. Finite numerical gate failures return the
diagnostic candidate, so a caller may retain rejected evidence without making
it a committable history.  A later root trial must use the same accepted parent.
"""
    validate_interval_energy_state(previous)
    row = _energy_record(energy)
    records = deepcopy(previous['records'])+[row]
    summary = _energy_summary(records)
    gates = _energy_gates(records, summary)
    accepted = all(gates.values())
    candidate = dict(contract=ENERGY_CONTRACT,
                     absolute_allowance_J=ENERGY_ABSOLUTE_J,
                     relative_allowance=ENERGY_RELATIVE,
                     records=records, **summary)
    return dict(accepted=accepted, gates=gates,
                reasons=[name for name, passed in gates.items() if not passed],
                step_limit_J=ENERGY_ABSOLUTE_J+ENERGY_RELATIVE*row['absolute_work_scale_J'],
                candidate_interval=candidate,
                next_state=deepcopy(candidate) if accepted else None)


def require_interval_energy(previous, energy):
    """Return a detached accepted accumulator or reject with ValueError."""
    result = assess_interval_energy(previous, energy)
    if not result['accepted']:
        raise ValueError('registered interval energy failed: '+', '.join(result['reasons']))
    return result['next_state']


def _pose(tool_pose):
    origin = np.asarray(tool_pose['origin'], dtype=float)
    rotation = np.asarray(tool_pose['rotation'], dtype=float)
    if (origin.shape != (3,) or rotation.shape != (3, 3) or
            not np.isfinite(origin).all() or not np.isfinite(rotation).all()):
        raise ValueError('invalid finite tool pose')
    if (not np.allclose(rotation.T@rotation, np.eye(3), rtol=0., atol=1e-12) or
            not np.allclose(rotation[:, 2], [0., 0., 1.], rtol=0., atol=1e-12) or
            abs(np.linalg.det(rotation)-1.) > 1e-12):
        raise ValueError('bore event observer requires proper rotation about z')
    return origin, rotation


def _particles(position, particle_ids, component_ids, masses):
    x = np.asarray(position, dtype=float)
    ids = np.asarray(particle_ids)
    components = np.asarray(component_ids)
    m = np.asarray(masses, dtype=float)
    if (x.ndim != 2 or x.shape[1:] != (3,) or ids.shape != (len(x),) or
            components.shape != ids.shape or m.shape != ids.shape or
            not np.isfinite(x).all() or not np.isfinite(m).all() or np.any(m <= 0)):
        raise ValueError('invalid finite particle geometry/mass')
    if (not np.issubdtype(ids.dtype, np.integer) or
            not np.issubdtype(components.dtype, np.integer) or
            len(np.unique(ids)) != len(ids)):
        raise ValueError('unique integer persistent IDs and integer components required')
    return x, ids, components, m


def initial_bore_crossing_state(position, particle_ids, component_ids, masses,
                                time_s, tool_pose, bore_width):
    """Freeze initially above-shoulder bar IDs; initial lower material is ineligible."""
    x, ids, components, m = _particles(position, particle_ids, component_ids, masses)
    t = _finite(time_s, 'initial observation time')
    width = _finite(bore_width, 'bore width')
    if width <= 0:
        raise ValueError('positive bore width required')
    origin, _ = _pose(tool_pose)
    ordering = np.argsort(ids)
    eligible = (components == 1) & (x[:, 2] > origin[2])
    return dict(contract=CROSSING_CONTRACT, initial_time_s=t, last_time_s=t,
                bore_width_m=width, particle_ids=ids[ordering].astype(int).tolist(),
                component_ids=components[ordering].astype(int).tolist(),
                masses_kg=m[ordering].tolist(),
                eligible_ids=sorted(ids[eligible].astype(int).tolist()),
                crossing_ids=[], crossing_mass_kg=0., events=[])


def validate_bore_crossing_state(state):
    required = {'contract', 'initial_time_s', 'last_time_s', 'bore_width_m',
                'particle_ids', 'component_ids', 'masses_kg', 'eligible_ids',
                'crossing_ids', 'crossing_mass_kg', 'events'}
    if set(state) != required or state['contract'] != CROSSING_CONTRACT:
        raise ValueError('bore event state fields or contract mismatch')
    ids, components = np.asarray(state['particle_ids']), np.asarray(state['component_ids'])
    _, ids, components, mass = _particles(np.zeros((len(ids), 3)), ids, components,
                                         state['masses_kg'])
    if not np.array_equal(ids, np.sort(ids)):
        raise ValueError('bore event IDs must be sorted')
    ti, t = (_finite(state[k], k) for k in ('initial_time_s', 'last_time_s'))
    width = _finite(state['bore_width_m'], 'bore width')
    if t < ti or width <= 0:
        raise ValueError('invalid bore event time or width')
    for name in ('eligible_ids', 'crossing_ids'):
        values = state[name]
        if values != sorted(set(values)) or any(type(v) is not int for v in values):
            raise ValueError('invalid '+name)
    if (not set(state['eligible_ids']).issubset(set(ids[components == 1])) or
            not set(state['crossing_ids']).issubset(set(state['eligible_ids']))):
        raise ValueError('bore event identity/eligibility mismatch')
    mass_by_id = dict(zip(ids.astype(int).tolist(), mass.tolist()))
    expected_mass = math.fsum(mass_by_id[i] for i in state['crossing_ids'])
    if state['crossing_mass_kg'] != expected_mass:
        raise ValueError('bore crossing mass mismatch')
    if not isinstance(state['events'], list) or len(state['events']) != len(state['crossing_ids']):
        raise ValueError('bore event count mismatch')
    event_ids = []
    for event in state['events']:
        if set(event) != {'particle_id', 'time_s', 'step_fraction', 'position_m',
                          'tool_local_position_m', 'mass_kg'}:
            raise ValueError('invalid bore event record fields')
        pid = event['particle_id']
        if type(pid) is not int or pid not in state['crossing_ids']:
            raise ValueError('invalid bore event particle ID')
        et = _finite(event['time_s'], 'event time')
        alpha = _finite(event['step_fraction'], 'event fraction')
        local = np.asarray(event['tool_local_position_m'], dtype=float)
        location = np.asarray(event['position_m'], dtype=float)
        if (et < ti or et > t or alpha < 0 or alpha > 1 or
                local.shape != (3,) or location.shape != (3,) or
                not np.isfinite(local).all() or not np.isfinite(location).all() or
                abs(local[2]) > 1e-12 or np.any(np.abs(local[:2]) >= width/2) or
                event['mass_kg'] != mass_by_id[pid]):
            raise ValueError('invalid through-square event evidence')
        event_ids.append(pid)
    if sorted(event_ids) != state['crossing_ids']:
        raise ValueError('bore events do not match unique crossing IDs')
    return True


def observe_square_bore_crossings(state, x0, x1, particle_ids, component_ids,
                                   masses, time0_s, time1_s, pose_at_time):
    """Observe the piecewise-linear accepted particle path, without moving it.

Tool rotation/translation is evaluated at the actual plane-intersection time,
not at either endpoint.  Shoulder z must be fixed or translate linearly on the
step; general tilted/accelerating planes are deliberately unsupported.  This
is centre-through-aperture evidence, not a deposition or contact certificate.
"""
    validate_bore_crossing_state(state)
    a, ids, components, m = _particles(x0, particle_ids, component_ids, masses)
    b, _, _, _ = _particles(x1, particle_ids, component_ids, masses)
    t0, t1 = _finite(time0_s, 'start time'), _finite(time1_s, 'end time')
    if t0 != state['last_time_s'] or t1 <= t0:
        raise ValueError('bore observation must continue accepted time')
    ordering = np.argsort(ids)
    if (ids[ordering].astype(int).tolist() != state['particle_ids'] or
            components[ordering].astype(int).tolist() != state['component_ids'] or
            m[ordering].tolist() != state['masses_kg']):
        raise ValueError('persistent IDs, components or mass changed')
    o0, _ = _pose(pose_at_time(t0))
    o1, _ = _pose(pose_at_time(t1))
    om, _ = _pose(pose_at_time(.5*(t0+t1)))
    if abs(om[2]-.5*(o0[2]+o1[2])) > 1e-14:
        raise ValueError('nonlinear shoulder z motion is not supported')
    d0, d1 = a[:, 2]-o0[2], b[:, 2]-o1[2]
    eligible = set(state['eligible_ids'])-set(state['crossing_ids'])
    candidates = np.flatnonzero((d0 >= 0) & (d1 < 0))
    events = []
    for index in candidates:
        pid = int(ids[index])
        if pid not in eligible:
            continue
        fraction = float(d0[index]/(d0[index]-d1[index]))
        et = t0+fraction*(t1-t0)
        location = a[index]+fraction*(b[index]-a[index])
        origin, rotation = _pose(pose_at_time(et))
        local = rotation.T@(location-origin)
        if abs(local[2]) > 1e-12:
            raise ValueError('interpolated point is not on actual shoulder plane')
        if np.all(np.abs(local[:2]) < .5*state['bore_width_m']):
            events.append(dict(particle_id=pid, time_s=float(et),
                               step_fraction=fraction, position_m=location.tolist(),
                               tool_local_position_m=local.tolist(), mass_kg=float(m[index])))
    events.sort(key=lambda e: (e['time_s'], e['particle_id']))
    result = deepcopy(state)
    result['last_time_s'] = t1
    result['events'].extend(events)
    result['crossing_ids'] = sorted(set(state['crossing_ids']) |
                                    {e['particle_id'] for e in events})
    result['crossing_mass_kg'] = math.fsum(e['mass_kg'] for e in result['events'])
    validate_bore_crossing_state(result)
    return dict(new_events=deepcopy(events), crossing_ids=result['crossing_ids'].copy(),
                crossing_mass_kg=result['crossing_mass_kg'],
                new_crossing_mass_kg=math.fsum(e['mass_kg'] for e in events),
                next_state=result)
