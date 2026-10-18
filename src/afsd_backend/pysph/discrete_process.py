"""DG-normal variant of the independently initialized hot quasi-static model.

The original endpoint material map and its raw energy defect are retained.
Only this variant accepts the new contact force and interval energy history.
"""
import numpy as np
from . import quasistatic_process as original
from .mechanics import positions
from .interval_acceptance import (initial_interval_energy_state,
    require_interval_energy, validate_interval_energy_state,
    initial_bore_crossing_state, observe_square_bore_crossings,
    validate_bore_crossing_state)

CLASSIFICATION = 'HOT_QUASISTATIC_CHAIN_RULE_NORMAL_CONTACT_APPROXIMATE'


def independent_initial_state(p):
    pa, state = original.independent_initial_state(p)
    state.update(classification=CLASSIFICATION,
                 interval_energy=initial_interval_energy_state())
    state['bore_observation'] = initial_bore_crossing_state(positions(pa), pa.particle_id,
        pa.component_id, pa.m, 0., original.pose(p, 0.), p['bore_width'])
    return pa, state


class DiscreteContactTrial(original.QuasistaticTrial):
    def __init__(self, pa, state, p):
        from .discrete_contact import ChainRuleDiscreteContact
        if state.get('classification') != CLASSIFICATION:
            raise ValueError('DG variant requires its own independent history')
        validate_interval_energy_state(state['interval_energy'])
        if state['interval_energy']['steps'] != state['step']:
            raise ValueError('interval history cannot be reset between accepted steps')
        validate_bore_crossing_state(state['bore_observation'])
        # Reuse the unchanged material, geometry, thermal and equilibrium map.
        # The temporary adapter is not published or used to restore old states.
        adapter = dict(state, classification=original.CLASSIFICATION)
        super().__init__(pa, adapter, p)
        self.state = state
        self.contact = ChainRuleDiscreteContact(tool_radius=p['tool_radius'],
            bore_width=p['bore_width'], tool_height=p['tool_top']-p['layer'],
            particle_radius=p['dx']/2, normal_stiffness=p['normal_stiffness'],
            tangential_stiffness=p['normal_stiffness'], friction=p['friction'],
            base_z=-p['substrate_depth'])


def accepted_copy(trial, result):
    """Recompute the supplied displacement, then enforce all-prefix energy gates."""
    fresh = trial.evaluate(result['displacement'], result['dt_s'])
    energy = original.energy_diagnostics(trial, fresh)
    next_interval = require_interval_energy(trial.state['interval_energy'], energy)
    pa, state = original.accepted_copy(trial, fresh)
    observation = observe_square_bore_crossings(trial.state['bore_observation'],
        trial.x0, positions(pa), pa.particle_id, pa.component_id, pa.m,
        trial.state['time_s'], state['time_s'], lambda t: original.pose(trial.p, t))
    state.update(interval_energy=next_interval, bore_observation=observation['next_state'],
        crossing_ids=observation['crossing_ids'], crossing_mass_kg=observation['crossing_mass_kg'])
    return pa, state, fresh
