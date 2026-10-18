"""Passive evidence of outlet passage and shoulder spreading, using arrays only.

No solver, particle update, force, heat, label or restart operation is imported.
Stage names are observational; even all stages do not certify deposition.
"""
import math
import numpy as np

CONTRACT = 'PASSIVE_ELIGIBLE_OUTLET_SHOULDER_STAGE_OBSERVER_V1'


def _array(value, shape, name, dtype=float):
    result = np.asarray(value, dtype=dtype)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError('invalid ' + name)
    return result


def _integer_ids(value, name):
    array = np.asarray(value)
    if array.ndim != 1 or array.dtype.kind not in 'iu' or len(np.unique(array)) != len(array):
        raise ValueError('unique integer ' + name + ' required')
    return array


def _longest_sampled_run(time, mask):
    longest = 0.; first = None
    for t, present in zip(time, mask):
        if present:
            if first is None:
                first = t
            longest = max(longest, float(t-first))
        else:
            first = None
    return longest


def observe_deposition_stages(time_s, particle_id, component_id, mass_kg,
        position_m, plastic_strain, tool_origin_m, tool_rotation,
        eligible_ids, crossing_events, *, dx_m, bore_width_m, tool_radius_m,
        shoulder_contact=None, shoulder_contact_available=None):
    """Return new JSON-compatible observations without altering any argument.

    Frame IDs must have one persistent ordering. ``eligible_ids`` and events
    come from the immutable initially-above-shoulder observer. Events cannot be
    inferred from an initially lower particle. If the input starts after an
    event, its first supplied frame is explicitly the plastic-strain baseline.
    Shoulder flags, if supplied, must be actual saved contact diagnostics.
    """
    ids = _integer_ids(particle_id, 'particle IDs'); n = len(ids)
    time = np.asarray(time_s, dtype=float); frames = len(time)
    if time.ndim != 1 or frames < 1 or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError('strictly increasing actual sample times required')
    component = np.asarray(component_id)
    if component.shape != (n,) or component.dtype.kind not in 'iu':
        raise ValueError('integer component IDs required')
    mass = _array(mass_kg, (n,), 'mass')
    if np.any(mass <= 0):
        raise ValueError('positive mass required')
    x = _array(position_m, (frames,n,3), 'positions')
    ep = _array(plastic_strain, (frames,n), 'plastic strain')
    origin = _array(tool_origin_m, (frames,3), 'tool origins')
    rotation = _array(tool_rotation, (frames,3,3), 'tool rotations')
    if (np.max(np.abs(rotation.swapaxes(-1,-2)@rotation-np.eye(3))) > 1e-12 or
            np.max(np.abs(np.linalg.det(rotation)-1)) > 1e-12):
        raise ValueError('proper orthogonal rotations required')
    dx, width, radius = float(dx_m), float(bore_width_m), float(tool_radius_m)
    if not np.isfinite([dx,width,radius]).all() or min(dx,width,radius) <= 0 or radius <= width/np.sqrt(2):
        raise ValueError('positive finite square-bore annular shoulder geometry required')
    eligible = _integer_ids(eligible_ids, 'eligible IDs')
    index = {int(pid): k for k,pid in enumerate(ids)}
    if any(int(pid) not in index or component[index[int(pid)]] != 1 for pid in eligible):
        raise ValueError('eligible IDs must be supplied bar particle IDs')
    eligible_set = set(int(pid) for pid in eligible)
    active = np.zeros((frames,n), dtype=bool)
    available = np.zeros(frames, dtype=bool)
    if shoulder_contact is not None:
        active = np.asarray(shoulder_contact)
        if active.shape != (frames,n) or active.dtype.kind != 'b':
            raise ValueError('actual boolean shoulder contact array required')
        available = np.ones(frames, dtype=bool)
        if shoulder_contact_available is not None:
            available = np.asarray(shoulder_contact_available)
            if available.shape != (frames,) or available.dtype.kind != 'b':
                raise ValueError('boolean contact availability required')
    elif shoulder_contact_available is not None:
        raise ValueError('contact availability requires actual contact array')

    local = np.einsum('fni,fij->fnj', x-origin[:,None,:], rotation)
    half = width/2; guard = dx/2
    square_offset = np.abs(local[:,:,:2])-half
    outside_square_distance = np.linalg.norm(np.maximum(square_offset,0),axis=2)
    signed_square_distance = np.where(np.max(square_offset,axis=2)>0,
        outside_square_distance, np.max(square_offset,axis=2))
    radial = np.linalg.norm(local[:,:,:2],axis=2)
    crossed = np.zeros((frames,n),dtype=bool)
    ep_delta = np.zeros((frames,n))
    event_by_id = {}; baseline_by_id = {}
    for event in crossing_events:
        pid = event['particle_id']
        if type(pid) is not int or pid not in eligible_set or pid in event_by_id:
            raise ValueError('crossing must be a unique originally eligible ID')
        et = float(event['time_s']); location = np.asarray(event['tool_local_position_m'],dtype=float)
        if (not np.isfinite(et) or et > time[-1] or location.shape != (3,) or
                not np.isfinite(location).all() or abs(location[2]) > 1e-12 or
                np.any(np.abs(location[:2]) >= half) or float(event['mass_kg']) != mass[index[pid]]):
            raise ValueError('invalid actual through-square event evidence')
        k = index[pid]; sampled = np.flatnonzero(time >= et)
        if not len(sampled):
            raise ValueError('event has no subsequent accepted sample')
        first = int(sampled[0]); crossed[first:,k] = True
        ep_delta[first:,k] = ep[first:,k]-ep[first,k]
        baseline_by_id[pid] = first; event_by_id[pid] = event
    below = local[:,:,2] < 0
    outside = signed_square_distance > 0
    under_shoulder = below & outside & (radial < radius)
    # Half one initial particle spacing is an explicit observation guard,
    # not a converged material-support radius or an added mechanical boundary.
    separated = (below & (signed_square_distance >= guard) & (radial <= radius-guard))
    ep_roundoff = 128*np.finfo(float).eps*np.maximum(1., np.abs(ep))
    plastic = ep_delta > ep_roundoff
    masks = dict(eligible_centre_crossed=crossed,
        crossed_currently_below_outlet=crossed & below,
        crossed_centre_under_annular_shoulder=crossed & under_shoulder,
        crossed_laterally_half_spacing_separated=crossed & separated,
        crossed_separated_with_post_crossing_plastic_increment=crossed & separated & plastic,
        crossed_actual_shoulder_contact=crossed & active & available[:,None],
        crossed_nominal_radius_fully_below_outlet=crossed & (local[:,:,2] <= -guard))
    history = []
    for frame,t in enumerate(time):
        row = dict(time_s=float(t),contact_observation_available=bool(available[frame]))
        for name, mask in masks.items():
            row[name+'_count'] = int(np.count_nonzero(mask[frame]))
            row[name+'_mass_kg'] = math.fsum(mass[mask[frame]].tolist())
        history.append(row)
    rows = []
    for pid in sorted(event_by_id):
        k = index[pid]; first = baseline_by_id[pid]
        row = dict(particle_id=pid,mass_kg=float(mass[k]),crossing_time_s=float(event_by_id[pid]['time_s']),
            post_crossing_baseline_sample_time_s=float(time[first]),
            baseline_is_first_available_accepted_sample_not_interpolated_event=True,
            position_m=x[-1,k].tolist(),tool_local_position_m=local[-1,k].tolist(),
            outlet_height_m=float(local[-1,k,2]),radial_distance_m=float(radial[-1,k]),
            signed_distance_from_square_projection_m=float(signed_square_distance[-1,k]),
            distance_to_square_projection_exit_m=float(max(-signed_square_distance[-1,k],0.)),
            distance_to_half_spacing_separation_m=float(max(guard-signed_square_distance[-1,k],0.)),
            shoulder_outer_margin_m=float(radius-radial[-1,k]),
            plastic_strain_at_baseline=float(ep[first,k]),plastic_strain_final=float(ep[-1,k]),
            post_crossing_plastic_increment=float(ep_delta[-1,k]),
            maximum_sampled_consecutive_separated_duration_s=_longest_sampled_run(time,masks['crossed_laterally_half_spacing_separated'][:,k]),
            actual_shoulder_contact_sample_count=int(np.count_nonzero(active[:,k] & available & crossed[:,k])))
        for name, mask in masks.items():
            row[name] = bool(mask[-1,k])
        rows.append(row)
    eligible_mask = np.isin(ids,eligible)
    initial_ineligible_bar = (component == 1) & ~eligible_mask
    final_crossed = crossed[-1]
    return dict(contract=CONTRACT,classification='PASSIVE_SAVED_ARRAY_DIAGNOSTIC_NOT_PROCESS_VALIDATION',
        physical_solver_evaluations=0,frame_count=frames,particle_count=n,
        time_start_s=float(time[0]),time_end_s=float(time[-1]),dx_m=dx,
        lateral_resolution_guard_m=guard,guard_basis='HALF_INITIAL_PARTICLE_SPACING_OBSERVATION_ONLY_NOT_SPATIAL_CONVERGENCE',
        eligible_initial_particle_count=len(eligible),eligible_initial_mass_kg=math.fsum(mass[eligible_mask].tolist()),
        eligible_remaining_uncrossed_count=int(np.count_nonzero(eligible_mask & ~final_crossed)),
        eligible_remaining_uncrossed_mass_kg=math.fsum(mass[eligible_mask & ~final_crossed].tolist()),
        initial_ineligible_bar_count=int(np.count_nonzero(initial_ineligible_bar)),
        initial_ineligible_bar_mass_kg=math.fsum(mass[initial_ineligible_bar].tolist()),
        ineligible_bar_currently_under_shoulder_count=int(np.count_nonzero(initial_ineligible_bar & under_shoulder[-1])),
        ineligible_bar_actual_shoulder_contact_count=int(np.count_nonzero(initial_ineligible_bar & active[-1] & available[-1])),
        final=history[-1],history=history,crossed_particle_rows=rows,
        complete_deposition_established=False,
        limitations=['Events and restart validity must be independently checked against the accepted chain.',
            'Consecutive stored samples do not prove continuous persistence between samples.',
            'A centre crossing, an initial lower particle or one shoulder contact cannot establish deposition.',
            'Half-spacing guards describe geometry at this resolution; they do not validate SPH spatial convergence.',
            'Plastic increments start at the first available accepted post-event sample; no event ep is invented.',
            'No bonding, steady deposition rate, substrate attachment, experimental or spatial convergence proof is provided.'])
