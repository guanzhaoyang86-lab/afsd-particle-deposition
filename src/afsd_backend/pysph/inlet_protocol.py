"""Unintegrated inlet inventory protocol, not a mass-flow boundary condition.

All operations are detached candidates. No particle is moved, no state is reset,
and no material/contact/thermal updater is called. A solver must separately
accept and atomically publish a new PA plus registry after rebuilding neighbours.
"""
import copy
import hashlib
import json

import numpy as np


CLASSIFICATION = 'INLET_PROTOCOL_CANDIDATE_NOT_ACCEPTED_NOT_INTEGRATED'
ACKNOWLEDGEMENT = 'NEW_REGISTRY_ONLY_NO_PRIOR_INLET'
REQUIRED = {
    'x': 'm', 'y': 'm', 'z': 'm', 'x0': 'm', 'y0': 'm', 'z0': 'm',
    'u': 'm/s', 'v': 'm/s', 'w': 'm/s', 'm': 'kg', 'rho': 'kg/m3',
    'h': 'm', 'temperature': 'K', 'S': 'Pa', 'ep': '1',
    'plastic_work': 'J', 'return_dissipation': 'J',
    'particle_id': '1', 'material_id': '1', 'component_id': '1',
    'tracer_type': '1', 'wholebar_tag': '1', 'actuator': '1',
}
DERIVED_STATUS = 'INITIAL_CONDITION'


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _sealed(value):
    result = copy.deepcopy(value)
    result.pop('registry_sha256', None)
    result['registry_sha256'] = hashlib.sha256(_json(result).encode()).hexdigest()
    return result


def _array_digest(pa):
    """Exact all-property hash aligned by persistent ID; sorting is not a reset."""
    n = pa.get_number_of_particles()
    if n <= 0 or len(pa.particle_id) != n:
        raise ValueError('nonempty real-particle inventory required')
    ids = np.asarray(pa.particle_id)
    if ids.dtype.kind not in 'iu' or np.any(ids <= 0) or len(np.unique(ids)) != n:
        raise ValueError('unique positive integer persistent IDs required')
    order = np.argsort(ids)
    digest = hashlib.sha256()
    digest.update(pa.name.encode())
    for key, prop in sorted(pa.properties.items()):
        value = prop.get_npy_array()
        stride = int(pa.stride.get(key, 1))
        if value.size != n*stride or not np.isfinite(value).all():
            raise ValueError('finite real-particle property layout required: '+key)
        digest.update(_json([key, value.dtype.str, stride]).encode())
        digest.update(value.reshape(n, stride)[order].tobytes())
    for key, value in sorted(pa.constants.items()):
        array = np.asarray(value)
        if not np.isfinite(array).all():
            raise ValueError('nonfinite constant: '+key)
        digest.update(_json([key, array.dtype.str, list(array.shape)]).encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


def initialize_registry(pa, provenance, *, acknowledgement):
    """Explicit one-time registry creation; never reconstruct this on restart."""
    if acknowledgement != ACKNOWLEDGEMENT or not isinstance(provenance, str) or not provenance.strip():
        raise ValueError('explicit new-registry acknowledgement and provenance required')
    digest = _array_digest(pa)
    highest = int(np.max(pa.particle_id))
    return _sealed(dict(version=1, classification=CLASSIFICATION, accepted=False,
        integration_status='NOT_INTEGRATED', provenance=provenance,
        initial_max_id=highest, high_watermark=highest, issued_particle_count=0,
        transaction_count=0, event_ids=[], last_event_time_s=None,
        thermal_ledger_convention=None,
        particle_state_sha256=digest,
        cumulative=dict(mass_kg=0., momentum_kg_m_s=[0., 0., 0.],
                        sensible_heat_J=0., kinetic_energy_J=0., volume_m3=0.)))


def restore_registry(serialized_json, pa):
    """Roundtrip registry and its PA together, never infer allocator from max(ID)."""
    registry = json.loads(serialized_json)
    if not isinstance(registry, dict) or _sealed(registry) != registry:
        raise ValueError('registry checksum mismatch')
    if (registry.get('version') != 1 or registry.get('classification') != CLASSIFICATION
            or registry.get('accepted') is not False or registry.get('integration_status') != 'NOT_INTEGRATED'):
        raise ValueError('unsupported registry classification')
    for key in ('initial_max_id', 'high_watermark', 'issued_particle_count', 'transaction_count'):
        if type(registry.get(key)) is not int or registry[key] < 0:
            raise ValueError('invalid allocator integer: '+key)
    if registry['high_watermark'] != registry['initial_max_id']+registry['issued_particle_count']:
        raise ValueError('allocator high watermark mismatch')
    events = registry['event_ids']
    if len(events) != registry['transaction_count'] or len(set(events)) != len(events):
        raise ValueError('registry event history mismatch')
    if int(np.max(pa.particle_id)) > registry['high_watermark']:
        raise ValueError('unregistered persistent ID')
    if registry['particle_state_sha256'] != _array_digest(pa):
        raise ValueError('registry particle history mismatch')
    return copy.deepcopy(registry)


def _parameter(spec, key, unit):
    row = spec[key]
    if (row.get('unit') != unit or row.get('status') != 'DECLARED_ASSUMPTION'
            or not isinstance(row.get('source'), str) or not row['source'].strip()):
        raise ValueError('explicit parameter unit/source/status required: '+key)
    value = float(row['value'])
    if not np.isfinite(value) or value <= 0:
        raise ValueError('positive finite parameter required: '+key)
    return value


def propose_insertion(pa, registry, incoming_pa, specification):
    """Return ``(candidate_pa, candidate_registry, source_event)`` atomically.

    Incoming fields already exist and all are declared. Persistent ID must be a
    zero placeholder. Scratch UNRESOLVED is recorded, never treated as solved.
    Supplying velocity is an inlet initial condition, not a trajectory law.
    """
    old = restore_registry(_json(registry), pa)
    spec = copy.deepcopy(specification)
    event_id = spec.get('event_id')
    if not isinstance(event_id, str) or not event_id.strip() or event_id in old['event_ids']:
        raise ValueError('nonempty unique inlet event_id required')
    if not isinstance(spec.get('source'), str) or not spec['source'].strip():
        raise ValueError('explicit incoming material source required')
    time = float(spec['time_s'])
    if (not np.isfinite(time) or time < 0 or
            (old['last_event_time_s'] is not None and time < old['last_event_time_s'])):
        raise ValueError('ordered finite inlet event time required')
    cp = _parameter(spec, 'heat_capacity', 'J/(kg K)')
    reference = _parameter(spec, 'reference_temperature', 'K')
    convention = dict(heat_capacity_J_kg_K=cp, reference_temperature_K=reference)
    if old['thermal_ledger_convention'] is not None and old['thermal_ledger_convention'] != convention:
        raise ValueError('cumulative thermal ledger convention cannot change')
    if pa.name != incoming_pa.name or set(pa.properties) != set(incoming_pa.properties):
        raise ValueError('incoming complete property schema/name must match')
    if not set(REQUIRED).issubset(pa.properties):
        raise ValueError('required physical/identity properties missing')
    declarations = spec.get('fields', {})
    if set(declarations) != set(pa.properties):
        raise ValueError('every incoming property requires explicit field declaration')
    n = incoming_pa.get_number_of_particles()
    if n <= 0:
        raise ValueError('nonempty incoming material required')
    unresolved = []
    for key, prop in pa.properties.items():
        values = incoming_pa.properties[key].get_npy_array()
        stride = int(pa.stride.get(key, 1))
        if (values.dtype != prop.get_npy_array().dtype or
                int(incoming_pa.stride.get(key, 1)) != stride or values.size != n*stride):
            raise ValueError('incoming schema dtype/stride/count mismatch: '+key)
        if not np.isfinite(values).all():
            raise ValueError('nonfinite incoming property: '+key)
        row = declarations[key]
        if not isinstance(row.get('source'), str) or not row['source'].strip() or not row.get('unit'):
            raise ValueError('field unit/source required: '+key)
        if key in REQUIRED and (row['unit'] != REQUIRED[key] or row.get('status') != DERIVED_STATUS):
            raise ValueError('required field unit/initial condition missing: '+key)
        if row.get('status') not in (DERIVED_STATUS, 'UNRESOLVED'):
            raise ValueError('unsupported field status: '+key)
        if row['status'] == 'UNRESOLVED':
            unresolved.append(key)
    if set(pa.constants) != set(incoming_pa.constants) or any(
            not np.array_equal(pa.constants[k], incoming_pa.constants[k]) for k in pa.constants):
        raise ValueError('incoming constants must match existing particle array')
    if np.any(incoming_pa.particle_id != 0):
        raise ValueError('incoming IDs must be unallocated zero placeholders')
    for key in ('m', 'rho', 'h', 'temperature'):
        if np.any(getattr(incoming_pa, key) <= 0):
            raise ValueError('positive incoming '+key+' required')
    if np.any(incoming_pa.ep < 0):
        raise ValueError('negative incoming plastic history')
    if np.any(incoming_pa.component_id != 1) or not np.isin(incoming_pa.material_id, pa.material_id).all():
        raise ValueError('incoming deformable feedstock and known material ID required')
    maximum = np.iinfo(pa.particle_id.dtype).max
    if old['high_watermark']+n > maximum:
        raise OverflowError('persistent ID integer storage exhausted')
    incoming = incoming_pa.extract_particles(np.arange(n, dtype=np.uint32))
    incoming.particle_id[:] = np.arange(old['high_watermark']+1, old['high_watermark']+n+1,
                                      dtype=pa.particle_id.dtype)
    candidate = pa.extract_particles(np.arange(pa.get_number_of_particles(), dtype=np.uint32))
    candidate.append_parray(incoming)
    candidate.set_output_arrays(list(candidate.properties))
    mass = np.asarray(incoming.m)
    velocity = np.column_stack((incoming.u, incoming.v, incoming.w))
    values = dict(mass_kg=float(np.sum(mass)),
        momentum_kg_m_s=np.sum(mass[:, None]*velocity, axis=0).tolist(),
        sensible_heat_J=float(np.sum(mass*cp*(incoming.temperature-reference))),
        kinetic_energy_J=float(np.sum(.5*mass*np.sum(velocity*velocity, axis=1))),
        volume_m3=float(np.sum(mass/incoming.rho)))
    event = dict(classification=CLASSIFICATION, accepted=False, integration_status='NOT_INTEGRATED',
        event_id=event_id, time_s=time, source=spec['source'],
        allocated_particle_ids=incoming.particle_id.tolist(), field_declarations=declarations,
        unresolved_properties=sorted(unresolved), source_increment=values,
        heat_capacity=spec['heat_capacity'], reference_temperature=spec['reference_temperature'],
        parent_particle_state_sha256=old['particle_state_sha256'],
        candidate_particle_state_sha256=_array_digest(candidate),
        pressure_flow_work='NOT_COMPUTED', elastic_energy='NOT_COMPUTED',
        mechanical_reaction_work='NOT_COMPUTED', total_energy_closure='NOT_CERTIFIED',
        transport_crossing='NOT_MEASURED_INSERTION_IS_NOT_CROSSING',
        required_before_physical_acceptance=['geometry_overlap_check', 'NNPS_rebuild',
            'boundary_contact_consistency', 'coupled_acceptance', 'atomic_PA_registry_publication'])
    next_registry = copy.deepcopy(old)
    next_registry.update(high_watermark=old['high_watermark']+n,
        issued_particle_count=old['issued_particle_count']+n,
        transaction_count=old['transaction_count']+1,
        event_ids=old['event_ids']+[event_id], last_event_time_s=time,
        thermal_ledger_convention=convention,
        particle_state_sha256=event['candidate_particle_state_sha256'])
    for key, value in values.items():
        total = np.asarray(old['cumulative'][key])+np.asarray(value)
        next_registry['cumulative'][key] = total.tolist()
    _json(event)  # Fail if even finite inputs overflow the explicit source ledger.
    return candidate, _sealed(next_registry), event
