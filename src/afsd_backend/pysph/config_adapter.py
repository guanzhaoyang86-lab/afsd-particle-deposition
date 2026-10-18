"""Fail-closed DEBUG adapter; paper configurations are never approximated implicitly."""
from afsd_scaffold.config import load_resolved_case, assert_physics_ready, ConfigurationError

UNITS = {"dx":"m","h_ratio":"1","rho0":"kg/m^3","shear":"Pa","poisson":"1","width":"m","bar_length":"m","substrate_width":"m","substrate_depth":"m","bore_width":"m","tool_radius":"m","tool_top":"m","temperature":"K","base_temperature":"K","tool_temperature":"K","conductivity":"W/(m*K)","cp":"J/(kg*K)","friction":"1","penalty_factor":"1","contact_conductance":"W/(m^2*K)","heat_partition":"1","plastic_heat_fraction":"1","a1":"Pa","a2":"Pa","a3":"1","b1":"1","b2":"1","c1":"1","c2":"1","room":"K","melt":"K","reference_rate":"1/s","rate_floor":"1/s","cfl":"1","density_limit":"1","penetration_limit":"1","energy_relative_limit":"1","energy_absolute_limit":"J","max_speed":"m/s","wall_seconds":"s","marker_low":"m","marker_high":"m","lift_speed":"m/s","dwell_duration":"s"}


UNITS.update({'actuator_depth': 'm', 'contact_reference_spacing': 'm', 'spatial_operator': '1'})


def load_settings(path):
    config = load_resolved_case(path)
    if config.get('intended_run_classification') not in ('DEBUG', 'APPROXIMATE'):
        assert_physics_ready(config)
        raise ConfigurationError('Only named DEBUG/APPROXIMATE host variants are implemented; PAPER disabled')
    required = ['debug.'+k for k in UNITS]
    required += ['process.'+k for k in ('tool_angular_speed', 'tool_traverse_speed', 'layer_height', 'feedstock_axial_speed')]
    assert_physics_ready(config, required, run_classification=config['intended_run_classification'])
    data = config['parameters']['debug']
    if set(data) != set(UNITS):
        raise ConfigurationError('unknown/missing DEBUG parameter')
    p = {}
    for k, unit in UNITS.items():
        if data[k]['unit'] != unit:
            raise ConfigurationError('wrong unit for '+k)
        p[k] = data[k]['value']
    for key, source in [('omega', 'tool_angular_speed'), ('traverse', 'tool_traverse_speed'), ('layer', 'layer_height'), ('feed', 'feedstock_axial_speed')]:
        expected_unit = {'omega': 'rad/s', 'traverse': 'm/s', 'layer': 'm', 'feed': 'm/s'}[key]
        if config['parameters']['process'][source]['unit'] != expected_unit:
            raise ConfigurationError('wrong process unit for '+source)
        p[key] = config['parameters']['process'][source]['value']
    if p['layer']/p['dx'] < 3-1e-12:
        raise ConfigurationError('DEBUG layer must have at least three particle spacings')
    if not (0 < p['friction'] < 1 and 0 <= p['heat_partition'] <= 1 and p['plastic_heat_fraction'] == 1):
        raise ConfigurationError('unsupported thermal/contact parameters')
    if not (p['room'] <= p['temperature'] < p['melt'] and p['dx'] > 0 and p['wall_seconds'] > 0):
        raise ConfigurationError('invalid DEBUG scale/domain')
    if p['dwell_duration'] != 0 or p['lift_speed'] != 0:
        raise ConfigurationError('hot-start probe has no dwell/lift implementation; nonzero request rejected')
    for key in ('rho0', 'shear', 'width', 'bar_length', 'substrate_width', 'substrate_depth', 'bore_width', 'tool_radius', 'tool_top', 'conductivity', 'cp', 'penalty_factor', 'contact_conductance', 'reference_rate', 'rate_floor', 'cfl', 'max_speed'):
        if p[key] <= 0: raise ConfigurationError('nonpositive '+key)
    if not 0 < p['poisson'] < .5: raise ConfigurationError('invalid Poisson ratio')
    if p['spatial_operator'] not in ('corrected_work_conjugate', 'paper_table1_raw_audit'):
        raise ConfigurationError('unknown spatial operator; no silent fallback')
    if not (0 < p['actuator_depth'] <= p['bar_length'] and p['contact_reference_spacing'] > 0):
        raise ConfigurationError('invalid physical actuator/contact reference length')
    return config, p, required
