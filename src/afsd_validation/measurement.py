"""Read-only proxy measurements under VALIDATION_CONTRACT v1.

This is deliberately NOT a cross-sectional oxide-area reconstruction. It
accepts the existing VTK reader's points/scalars, selects one passive label
scope, and reports clearly named count/volume proxies only. All lengths are
metres in LAB coordinates; all masses kg and densities kg/m^3. Unknown paper
ROI/section/start-stop intervals must be supplied as explicit assumptions or
the call fails. No experimental reference values are present in this module.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
import hashlib
import json
import math


class MeasurementError(ValueError):
    """Ambiguous/invalid sampling is not a zero measurement."""


_SPEC_KEYS = {
    'frame', 'label_scope', 'tracer_value', 'material_id', 'roi_x_m',
    'section_origin_m', 'section_normal', 'section_thickness_m',
    'startup_x_m', 'termination_x_m', 'region', 'roi_source', 'roi_status',
    'time_s', 'omega_rad_s', 'traverse_m_s', 'tool_axis_y_m', 'side_band_m',
}
_SCOPES = {'LOCAL': 'tracer_type', 'WHOLEBAR_DIAGNOSTIC': 'wholebar_tag'}


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise MeasurementError(f'{name} requires an explicit finite SI number')
    return float(value)


def _vector(value, length, name):
    if not isinstance(value, (list, tuple)) or len(value) != length:
        raise MeasurementError(f'{name} requires {length} explicit values')
    return tuple(_number(x, name) for x in value)


def _interval(value, name):
    a, b = _vector(value, 2, name)
    if not a < b:
        raise MeasurementError(f'{name} requires ordered, nonempty endpoints')
    return a, b


def _overlap(a, b):
    return max(a[0], b[0]) < min(a[1], b[1])


def validate_spec(spec: Mapping):
    """Strict protocol: missing values and alternative/mixed labels fail closed."""
    if not isinstance(spec, Mapping) or set(spec) != _SPEC_KEYS:
        missing = sorted(_SPEC_KEYS-set(spec)) if isinstance(spec, Mapping) else sorted(_SPEC_KEYS)
        extra = sorted(set(spec)-_SPEC_KEYS) if isinstance(spec, Mapping) else []
        raise MeasurementError(f'explicit sampling specification required; missing={missing}, unexpected={extra}')
    if spec['frame'] != 'LAB':
        raise MeasurementError('positions/ROI must be LAB; no silent rotating-frame transformation')
    if not isinstance(spec['label_scope'], str) or spec['label_scope'] not in _SCOPES:
        raise MeasurementError('select exactly one LOCAL or WHOLEBAR_DIAGNOSTIC label scope')
    if type(spec['tracer_value']) is not int or spec['tracer_value'] not in (1, 2):
        raise MeasurementError('select exactly one outer=1 or core=2 tracer value')
    if type(spec['material_id']) is not int or spec['material_id'] != 1:
        raise MeasurementError('this protocol samples AA6061-T6 material_id=1 only')
    roi = _interval(spec['roi_x_m'], 'roi_x_m')
    startup = _interval(spec['startup_x_m'], 'startup_x_m')
    termination = _interval(spec['termination_x_m'], 'termination_x_m')
    if _overlap(startup, termination):
        raise MeasurementError('startup and termination windows may not overlap')
    region = spec['region']
    if region == 'STEADY_ASSUMED':
        if _overlap(roi, startup) or _overlap(roi, termination):
            raise MeasurementError('steady-assumed ROI must exclude declared startup/termination intervals')
    elif region in ('STARTUP', 'TERMINATION'):
        window = startup if region == 'STARTUP' else termination
        if not (window[0] <= roi[0] < roi[1] <= window[1]):
            raise MeasurementError('ROI is outside its declared start/stop region')
    else:
        raise MeasurementError('region must explicitly be STARTUP, STEADY_ASSUMED or TERMINATION')
    if spec['roi_status'] not in ('ASSUMED', 'DERIVED', 'REFERENCE_CONFIRMED') or not isinstance(spec['roi_source'], str) or not spec['roi_source'].strip():
        raise MeasurementError('explicit ROI source and epistemic status required')
    _vector(spec['section_origin_m'], 3, 'section_origin_m')
    normal = _vector(spec['section_normal'], 3, 'section_normal')
    if abs(sum(x*x for x in normal)-1) > 1e-12:
        raise MeasurementError('section_normal must be an explicit unit LAB vector')
    if _number(spec['section_thickness_m'], 'section_thickness_m') <= 0:
        raise MeasurementError('positive section thickness required')
    if _number(spec['time_s'], 'time_s') < 0:
        raise MeasurementError('negative observation time')
    for name in ('omega_rad_s', 'traverse_m_s'):
        if _number(spec[name], name) == 0:
            raise MeasurementError('AS/RS undefined at zero rotation or traverse')
    if _number(spec['tool_axis_y_m'], 'tool_axis_y_m') != 0:
        raise MeasurementError('v1 freezes tool-axis y=0; no outcome-dependent recentering')
    if _number(spec['side_band_m'], 'side_band_m') != 0.01588:
        raise MeasurementError('v1 freezes each side analysis band at 15.88 mm')


def _ratio(numerator, denominator):
    return {
        'value': numerator/denominator if denominator != 0 else None,
        'status': 'DEFINED' if denominator != 0 else 'UNDEFINED',
        'numerator_AS': numerator, 'denominator_RS': denominator,
        'reason': None if denominator != 0 else 'zero RS denominator; no epsilon or experimental-zero inference',
    }


def measure_proxies(points: Sequence, scalar_fields: Mapping, spec: Mapping):
    """Measure R_count and R_volume in an explicitly supplied geometric slab.

    ROI x endpoints are half-open [x0,x1); slab and side-band edges are closed.
    y=0 belongs to neither AS nor RS and is reported separately. Particle cells
    intersecting the slab are NOT reconstructed: only their centres are selected.
    Thus even volume/thickness is not returned as an experimental area.
    """
    validate_spec(spec)
    label_field = _SCOPES[spec['label_scope']]
    required = ('particle_id', 'material_id', label_field, 'mass_kg', 'density_kg_per_m3')
    count = len(points)
    for field in required:
        if field not in scalar_fields or len(scalar_fields[field]) != count:
            raise MeasurementError(f'missing/misaligned field {field}; no zero fill')
    ids = scalar_fields['particle_id']
    if any(type(pid) is not int or pid < 0 for pid in ids) or len(set(ids)) != count:
        raise MeasurementError('particle IDs must be unique nonnegative integers')
    origin = spec['section_origin_m']; normal = spec['section_normal']
    lo, hi = spec['roi_x_m']; half = spec['section_thickness_m']/2
    side_ids = {'AS': [], 'RS': [], 'CENTERLINE': []}
    side_volumes = {'AS': [], 'RS': []}
    for i, raw in enumerate(points):
        xyz = _vector(raw, 3, f'position[{i}]')
        mass = _number(scalar_fields['mass_kg'][i], 'mass_kg')
        density = _number(scalar_fields['density_kg_per_m3'][i], 'density_kg_per_m3')
        if mass <= 0 or density <= 0:
            raise MeasurementError('positive mass and density required for physical volume')
        label = scalar_fields[label_field][i]
        material = scalar_fields['material_id'][i]
        if type(label) is not int or type(material) is not int:
            raise MeasurementError('integer passive label and material fields required')
        if material != spec['material_id'] or label != spec['tracer_value']:
            continue
        distance = sum((xyz[k]-origin[k])*normal[k] for k in range(3))
        if not (lo <= xyz[0] < hi and abs(distance) <= half and abs(xyz[1]) <= spec['side_band_m']):
            continue
        product = -spec['omega_rad_s']*xyz[1]*spec['traverse_m_s']
        side = 'AS' if product > 0 else 'RS' if product < 0 else 'CENTERLINE'
        side_ids[side].append(ids[i])
        if side != 'CENTERLINE':
            side_volumes[side].append(mass/density)
    frozen_spec = deepcopy(dict(spec))
    signature = hashlib.sha256(json.dumps(frozen_spec, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest()
    return {
        'classification': 'APPROXIMATE', 'contract_version': 'VALIDATION_CONTRACT_v1',
        'operator': 'CENTRE_SELECTED_SLAB_COUNT_AND_VOLUME_PROXIES_ONLY',
        'specification': frozen_spec, 'specification_sha256': signature,
        'label_field': label_field, 'steady_state_verified': False,
        'R_count': _ratio(len(side_ids['AS']), len(side_ids['RS'])),
        'R_volume': _ratio(math.fsum(side_volumes['AS']), math.fsum(side_volumes['RS'])),
        'volume_unit': 'm^3', 'ratio_unit': '1',
        'selected_particle_ids': {side: sorted(value) for side, value in side_ids.items()},
        'experimental_area_ratio': {
            'value': None, 'status': 'BLOCKED',
            'reason': 'no cross-sectional marker-area reconstruction or validated experimental sampling; R_count/R_volume are not Fig.8 area ratios',
        },
    }
