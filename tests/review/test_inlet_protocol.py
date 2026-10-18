"""Manufactured inlet bookkeeping references, never AFSD or transport tests."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from pysph.solver.utils import load

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import inlet_protocol as inlet
from afsd_backend.pysph.output import checkpoint
from afsd_backend.pysph.particles import make_array


def fixture():
    p = dict(dx=.1, h_ratio=1.3, rho0=2., width=1., marker_low=0., marker_high=1.,
             bar_length=2., actuator_depth=.2, temperature=300., traverse=.2, feed=.1)
    old = make_array([[0., 0., .1], [.1, .1, .2]], np.array([0, 1]), p)
    old.particle_id[:] = [7, 3]
    old.ep[:] = [.2, .3]; old.plastic_work[:] = [11., 17.]
    old.return_dissipation[:] = [5., 6.]
    old.S[:] = np.arange(18)+.5
    incoming = make_array([[.1, .2, .5], [-.1, 0., .7]], np.ones(2, int), p)
    incoming.particle_id[:] = 0
    incoming.m[:] = [2., 3.]; incoming.rho[:] = [2., 3.]
    incoming.u[:], incoming.v[:], incoming.w[:] = [1., -1.], [2., 0.], [0., 2.]
    incoming.temperature[:] = [310., 280.]
    incoming.ep[:] = [.04, .09]; incoming.plastic_work[:] = [2., 9.]
    incoming.return_dissipation[:] = [1., 3.]
    for pa in (old, incoming):
        pa.add_property('as_rs_tag', type='long', data=[-1, 1])
        pa.add_property('extra_history', stride=2, data=[1., 3., 7., 9.])
        pa.set_output_arrays(list(pa.properties))
    fields = {key: dict(unit=inlet.REQUIRED.get(key, '1'),
        source='MANUFACTURED_INITIAL_STATE_NOT_AA6061',
        status='INITIAL_CONDITION' if key in inlet.REQUIRED or key in ('as_rs_tag', 'extra_history')
               else 'UNRESOLVED') for key in incoming.properties}
    spec = dict(event_id='feed-1', time_s=1., source='MANUFACTURED_EXTERNAL_RESERVOIR', fields=fields,
        heat_capacity=dict(value=10., unit='J/(kg K)', source='MANUFACTURED_REFERENCE', status='DECLARED_ASSUMPTION'),
        reference_temperature=dict(value=300., unit='K', source='MANUFACTURED_REFERENCE', status='DECLARED_ASSUMPTION'))
    state = inlet.initialize_registry(old, 'MANUFACTURED_NEW_INVENTORY', acknowledgement=inlet.ACKNOWLEDGEMENT)
    return old, state, incoming, spec


def properties(pa):
    return {k: v.get_npy_array().copy() for k, v in pa.properties.items()}


class InletProtocolTests(unittest.TestCase):
    def test_independent_mass_momentum_sensible_heat_kinetic_and_volume_reference(self):
        old, state, inc, spec = fixture()
        candidate, updated, event = inlet.propose_insertion(old, state, inc, spec)
        # Written independently, not evaluated with the inlet implementation.
        reference = dict(mass_kg=5., momentum_kg_m_s=[-1., 4., 6.],
                         sensible_heat_J=-400., kinetic_energy_J=12.5, volume_m3=2.)
        for key, answer in reference.items():
            np.testing.assert_allclose(event['source_increment'][key], answer, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(updated['cumulative'][key], answer, rtol=1e-12, atol=1e-12)
        self.assertEqual(candidate.get_number_of_particles(), 4)
        self.assertLess(event['source_increment']['sensible_heat_J'], 0.)
        self.assertEqual(event['transport_crossing'], 'NOT_MEASURED_INSERTION_IS_NOT_CROSSING')
        self.assertFalse(event['accepted'])
        self.assertEqual(event['total_energy_closure'], 'NOT_CERTIFIED')

    def test_arrays_prior_history_incoming_and_metadata_do_not_mutate_or_alias(self):
        old, state, inc, spec = fixture()
        saved_old, saved_inc, saved_state, saved_spec = properties(old), properties(inc), copy.deepcopy(state), copy.deepcopy(spec)
        candidate, updated, event = inlet.propose_insertion(old, state, inc, spec)
        for key, value in saved_old.items():
            stride = old.stride.get(key, 1)
            np.testing.assert_array_equal(candidate.properties[key].get_npy_array()[:2*stride], value)
            np.testing.assert_array_equal(old.properties[key].get_npy_array(), value)
            np.testing.assert_array_equal(inc.properties[key].get_npy_array(), saved_inc[key])
            if key != 'particle_id':
                np.testing.assert_array_equal(candidate.properties[key].get_npy_array()[2*stride:], saved_inc[key])
        candidate.S[0] = 999.; candidate.x[0] = 999.
        updated['event_ids'].append('mutation'); event['field_declarations']['x']['source'] = 'mutated'
        np.testing.assert_array_equal(old.S, saved_old['S'])
        self.assertEqual(state, saved_state); self.assertEqual(spec, saved_spec)

    def test_monotonic_ids_and_json_restart_allocate_only_new_ids(self):
        old, state, inc, spec = fixture()
        a, state, event = inlet.propose_insertion(old, state, inc, spec)
        self.assertEqual(event['allocated_particle_ids'], [8, 9])
        state = inlet.restore_registry(json.dumps(state, allow_nan=False), a)
        spec.update(event_id='feed-2', time_s=2.)
        b, state, event = inlet.propose_insertion(a, state, inc, spec)
        self.assertEqual(event['allocated_particle_ids'], [10, 11])
        self.assertEqual(state['high_watermark'], 11)
        self.assertEqual(state['issued_particle_count'], 4)
        self.assertEqual(state['cumulative']['mass_kg'], 10.)
        self.assertEqual(state['cumulative']['sensible_heat_J'], -800.)
        self.assertEqual(len(np.unique(b.particle_id)), 6)

    def test_registered_arrays_may_reorder_without_identity_corruption(self):
        old, state, inc, spec = fixture()
        old = old.extract_particles(np.array([1, 0], dtype=np.uint32))
        state = inlet.restore_registry(json.dumps(state), old)
        candidate, _, _ = inlet.propose_insertion(old, state, inc, spec)
        self.assertEqual(candidate.particle_id.tolist(), [3, 7, 8, 9])
        self.assertEqual(candidate.ep.tolist(), [.3, .2, .04, .09])

    def test_checkpoint_pa_registry_roundtrip_and_continue_ids(self):
        old, state, inc, spec = fixture()
        candidate, state, _ = inlet.propose_insertion(old, state, inc, spec)
        with tempfile.TemporaryDirectory(prefix='afsd_inlet_protocol_') as tmp:
            path = Path(tmp)/'protocol_not_physics.npz'
            checkpoint(path, candidate, {'classification': inlet.CLASSIFICATION, 'registry': state,
                                        'accepted': False, 'restart_allowed': False})
            data = load(str(path)); restored = data['arrays']['solid']
            state = inlet.restore_registry(json.dumps(data['solver_data']['registry']), restored)
            for key in candidate.properties:
                np.testing.assert_array_equal(restored.properties[key].get_npy_array(), candidate.properties[key].get_npy_array())
            spec.update(event_id='feed-2', time_s=2.)
            _, _, event = inlet.propose_insertion(restored, state, inc, spec)
            self.assertEqual(event['allocated_particle_ids'], [10, 11])

    def test_tag_changes_never_change_injection_ledger(self):
        old, state, inc, spec = fixture()
        a, _, first = inlet.propose_insertion(old, state, inc, spec)
        inc.tracer_type[:] = [88, -4]; inc.wholebar_tag[:] = [6, 8]; inc.as_rs_tag[:] *= -1
        b, _, changed = inlet.propose_insertion(old, state, inc, spec)
        self.assertEqual(first['source_increment'], changed['source_increment'])
        for key in ('m', 'rho', 'temperature', 'u', 'v', 'w', 'S', 'ep', 'x', 'y', 'z'):
            np.testing.assert_array_equal(getattr(a, key), getattr(b, key))

    def test_source_ledger_equals_direct_inventory_differences(self):
        old, state, inc, spec = fixture()
        new, _, event = inlet.propose_insertion(old, state, inc, spec)
        # Direct differences are independent of the implementation's incoming-only sum.
        def totals(pa):
            return np.r_[pa.m.sum(), (pa.m*pa.u).sum(), (pa.m*pa.v).sum(), (pa.m*pa.w).sum(),
                         (10.*pa.m*(pa.temperature-300.)).sum(),
                         (.5*pa.m*(pa.u**2+pa.v**2+pa.w**2)).sum()]
        source = event['source_increment']
        expected = np.r_[source['mass_kg'], source['momentum_kg_m_s'], source['sensible_heat_J'], source['kinetic_energy_J']]
        np.testing.assert_allclose(totals(new)-totals(old), expected, rtol=1e-12, atol=1e-12)

    def test_duplicate_event_stale_registry_and_state_tampering_fail(self):
        old, state, inc, spec = fixture()
        new, latest, _ = inlet.propose_insertion(old, state, inc, spec)
        with self.assertRaisesRegex(ValueError, 'unique inlet event'): inlet.propose_insertion(new, latest, inc, spec)
        with self.assertRaises(ValueError): inlet.restore_registry(json.dumps(state), new)
        corrupted = copy.deepcopy(latest); corrupted['high_watermark'] = 7
        with self.assertRaisesRegex(ValueError, 'checksum'): inlet.restore_registry(json.dumps(corrupted), new)
        new.ep[0] += .01
        with self.assertRaisesRegex(ValueError, 'particle history'): inlet.restore_registry(json.dumps(latest), new)

    def test_removed_particle_is_not_a_license_to_reuse_id(self):
        old, state, inc, spec = fixture()
        new, latest, _ = inlet.propose_insertion(old, state, inc, spec)
        removed = new.extract_particles(np.array([0, 1, 2], dtype=np.uint32))
        with self.assertRaisesRegex(ValueError, 'particle history'): inlet.restore_registry(json.dumps(latest), removed)
        with self.assertRaises(ValueError): inlet.initialize_registry(removed, 'restart', acknowledgement='')

    def test_missing_or_unresolved_required_field_unit_and_source_fail(self):
        for mode in ('missing', 'unresolved', 'wrong-unit', 'no-source'):
            old, state, inc, spec = fixture()
            if mode == 'missing': del spec['fields']['S']
            elif mode == 'unresolved': spec['fields']['S']['status'] = 'UNRESOLVED'
            elif mode == 'wrong-unit': spec['fields']['S']['unit'] = 'MPa'
            else: spec['fields']['S']['source'] = ''
            with self.subTest(mode=mode), self.assertRaises(ValueError): inlet.propose_insertion(old, state, inc, spec)

    def test_scratch_fields_stay_explicitly_unsolved(self):
        old, state, inc, spec = fixture()
        inc.e[:] = [123., 456.]
        candidate, _, event = inlet.propose_insertion(old, state, inc, spec)
        self.assertIn('e', event['unresolved_properties'])
        np.testing.assert_array_equal(candidate.e[-2:], [123., 456.])
        self.assertEqual(event['source_increment']['sensible_heat_J'], -400.)
        self.assertEqual(event['elastic_energy'], 'NOT_COMPUTED')

    def test_invalid_particle_values_fail_without_changing_existing_state(self):
        for key, value in (('m', 0.), ('rho', -1.), ('h', 0.), ('temperature', 0.),
                           ('S', np.nan), ('ep', -.1), ('particle_id', 8), ('component_id', 0), ('material_id', 66)):
            old, state, inc, spec = fixture(); snapshot = properties(old)
            getattr(inc, key)[0] = value
            with self.subTest(key=key), self.assertRaises(ValueError): inlet.propose_insertion(old, state, inc, spec)
            for field, values in snapshot.items(): np.testing.assert_array_equal(old.properties[field].get_npy_array(), values)

    def test_extra_missing_property_schema_and_time_fail(self):
        old, state, inc, spec = fixture(); inc.add_property('unannounced')
        with self.assertRaisesRegex(ValueError, 'schema'): inlet.propose_insertion(old, state, inc, spec)
        old, state, inc, spec = fixture(); spec['time_s'] = np.nan
        with self.assertRaisesRegex(ValueError, 'event time'): inlet.propose_insertion(old, state, inc, spec)
        old, state, inc, spec = fixture(); spec['heat_capacity']['unit'] = 'kJ/(kg K)'
        with self.assertRaisesRegex(ValueError, 'parameter unit'): inlet.propose_insertion(old, state, inc, spec)

    def test_cumulative_heat_cannot_mix_reference_temperatures_or_material_cp(self):
        old, state, inc, spec = fixture()
        new, state, _ = inlet.propose_insertion(old, state, inc, spec)
        for key, value in (('heat_capacity', 12.), ('reference_temperature', 273.15)):
            different = copy.deepcopy(spec)
            different.update(event_id='feed-2', time_s=2.)
            different[key]['value'] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'thermal ledger convention'):
                inlet.propose_insertion(new, state, inc, different)


if __name__ == '__main__': unittest.main(verbosity=2)
