"""Synthetic receipt integrity tests; no physical solver or GPU evaluations."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import afsd_event_receipt_20261010 as receipt
from afsd_event_blocked_path_20261010 import audit_blocked_prefix


class FixtureTransaction:
    def __init__(self, root):
        self.root = root
        self.states = {}
        self.particles = {}

    def sha(self, path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    def read(self, path):
        return json.loads(Path(path).read_text(encoding='utf-8'))

    def load(self, path):
        digest=self.sha(path)
        return copy.deepcopy(self.particles[digest]), copy.deepcopy(self.states[digest])


class EventReceiptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.tx = FixtureTransaction(self.root)
        self.folder = self.root/'stage'
        self.folder.mkdir()
        self.meta = {'artifacts_sha256': {}, 'source_sha256': {'fixture_source': 'a'*64}}
        parent = dict(time_s=.2, step=2, diagnostic_only=True, accepted=False, restart_allowed=False,
            machine_history=[{'t0_s': 0., 't1_s': .1, 'dt_s': .1}, {'t0_s': .1, 't1_s': .2, 'dt_s': .1}],
            interval_energy={'records': [{'energy': 0}, {'energy': 1}]})
        parent['implementation_source_sha256']=copy.deepcopy(self.meta['source_sha256'])
        parent['parameters']=dict(tool_radius=2.,bore_width=1.,tool_top=1.,layer=0.)
        for motion in parent['machine_history']:
            motion.update(origin0_m=[0.,0.,0.],origin1_m=[0.,0.,0.],rotation0=np.eye(3).tolist(),rotation1=np.eye(3).tolist())
        self.parent = self.root/'parent.npz'
        self.checkpoint(self.parent, parent, artifact=False)
        self.meta.update(parent_checkpoint=self.rel(self.parent), parent_sha256=self.tx.sha(self.parent),
            parent_fingerprint=self.fingerprint({}, parent), parent_time_s=.2, dt_s=.2,
            screened_continuation_authorized=True,
            event_substeps={'schema': receipt.SCHEMA}, screened_history=[],event_protocol='SYNTHETIC_EVENT_PROTOCOL',
            event_tolerance_scale=1.,event_records={role:[] for role in receipt.ROLES})
        states = {}
        for role, dts in [('coarse', [.1, .1]), ('fine_1', [.025, .075]), ('fine_2', [1e-9, .1-1e-9])]:
            if role == 'fine_2':
                previous = copy.deepcopy(states['fine_1'])
                previous_path = self.folder/'fine_1/NOT_RESTARTABLE_candidate.npz'
            else:
                previous, previous_path = copy.deepcopy(parent), self.parent
            rows = []
            for i, dt in enumerate(dts):
                leaf = self.folder/role/f'leaf_{i+1:04d}'
                checkpoint = leaf/'NOT_RESTARTABLE_candidate.npz'
                current = copy.deepcopy(previous)
                current['time_s'] += dt
                current['step'] += 1
                current['machine_history'].append({'t0_s': previous['time_s'], 't1_s': current['time_s'], 'dt_s': dt,
                    'origin0_m':[0.,0.,0.],'origin1_m':[0.,0.,0.],
                    'rotation0':np.eye(3).tolist(),'rotation1':np.eye(3).tolist()})
                current['interval_energy']['records'].append({'energy': role+str(i)})
                self.checkpoint(checkpoint, current)
                observer = leaf/'observer.json'
                qualification = leaf/'candidate_qualification.json'
                self.write(observer, {'time_s': current['time_s']})
                self.write(qualification, {'status': receipt.QUALIFIED, 'accepted': False, 'restart_allowed': False})
                rows.append(dict(checkpoint=self.rel(checkpoint), sha256=self.tx.sha(checkpoint),
                    parent_checkpoint=self.rel(previous_path), parent_sha256=self.tx.sha(previous_path),
                    parent_time_s=previous['time_s'], parent_step=previous['step'], time_s=current['time_s'],
                    step=current['step'], dt_s=dt, observer_path=self.rel(observer), qualification_path=self.rel(qualification)))
                previous, previous_path = current, checkpoint
            alias = self.folder/role/'NOT_RESTARTABLE_candidate.npz'
            self.write(alias, previous_path.read_bytes())
            observer_alias = self.folder/role/'observer.json'
            self.write(observer_alias, (self.root/rows[-1]['observer_path']).read_bytes())
            self.meta['event_substeps'][role] = rows
            states[role] = previous
            if role != 'coarse':
                self.meta['screened_history'].append(dict(checkpoint=self.rel(alias), sha256=self.tx.sha(alias),
                    time_s=previous['time_s'], observer_path=self.rel(observer_alias)))
        fine = self.meta['screened_history']
        self.meta.update(screened_continuation_checkpoint=fine[-1]['checkpoint'],
            screened_continuation_sha256=fine[-1]['sha256'], screened_end_time_s=states['fine_2']['time_s'],
            screened_end_step=states['fine_2']['step'],
            pair_screen=dict(passed=True, global_screen={'pass': True}, local_screen={'passed': True},
                crossing_ids_match=True, fine_checkpoint_sha256=[r['sha256'] for r in fine],
                coarse_checkpoint_sha256=self.tx.sha(self.folder/'coarse/NOT_RESTARTABLE_candidate.npz')))
        self.write(self.folder/'pair_screen.json', self.meta['pair_screen'])

    def rel(self, path):
        return Path(path).relative_to(self.root).as_posix()

    def write(self, path, value, artifact=True):
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(value, bytes):
            path.write_bytes(value)
        else:
            path.write_text(json.dumps(value), encoding='utf-8')
        if artifact:
            self.meta['artifacts_sha256'][self.rel(path)] = self.tx.sha(path)

    def checkpoint(self, path, state, artifact=True):
        self.write(path, json.dumps({'state': state, 'origin': self.rel(path)}).encode(), artifact)
        self.tx.states[self.tx.sha(path)] = copy.deepcopy(state)
        self.tx.particles[self.tx.sha(path)] = self.particles(.995)

    @staticmethod
    def particles(z):
        return SimpleNamespace(x=np.array([.4,.6]),y=np.zeros(2),z=np.array([z,z]),particle_id=np.array([10,20]))

    @staticmethod
    def fingerprint(pa, state):
        return hashlib.sha256(json.dumps(state, sort_keys=True).encode()).hexdigest()

    @staticmethod
    def state_validator(state):
        if not (len(state['machine_history']) == state['step'] == len(state['interval_energy']['records'])):
            raise ValueError('synthetic full history mismatch')

    def validate(self):
        return receipt.validate_gpu_screened_receipt(self.tx, self.meta, self.folder,
            legacy_validator=lambda *args: self.fail('event receipt used old +2 validator'),
            state_validator=self.state_validator, parent_fingerprint=self.fingerprint)

    def row(self, role='fine_1', index=0):
        return self.meta['event_substeps'][role][index]

    def test_multiple_real_steps_and_tiny_positive_event_interval_pass(self):
        result = self.validate()
        self.assertEqual(result['fine_physical_steps'], 4)
        self.assertEqual(result['substep_counts'], {'coarse': 2, 'fine_1': 2, 'fine_2': 2})

    def test_ordinary_receipt_calls_frozen_validator_without_event_schema(self):
        calls = []
        ordinary = {'no_event_schema': True}
        result = receipt.validate_gpu_screened_receipt(self.tx, ordinary, self.folder,
            legacy_validator=lambda *args: calls.append(args) or 'frozen-result')
        self.assertEqual(result, 'frozen-result')
        self.assertEqual(calls, [(self.tx, ordinary, self.folder)])

    def test_unknown_empty_or_missing_role_schema_rejected(self):
        for mutation in ({'schema': 'unknown'}, {'schema': receipt.SCHEMA}, None):
            with self.subTest(schema=mutation):
                self.meta['event_substeps'] = mutation
                with self.assertRaises((ValueError, AssertionError)):
                    self.validate()

    def test_incomplete_leaf_row_rejected(self):
        del self.row()['qualification_path']
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_missing_or_changed_qualification_artifact_rejected(self):
        path = self.root/self.row()['qualification_path']
        path.write_text('{}', encoding='utf-8')
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_unqualified_leaf_cannot_continue(self):
        path = self.root/self.row()['qualification_path']
        self.write(path, {'status': 'FAIL', 'accepted': False, 'restart_allowed': False})
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_foreign_qualification_source_rejected(self):
        path = self.root/self.row()['qualification_path']
        self.write(path, {'status': receipt.QUALIFIED, 'accepted': False, 'restart_allowed': False,
                          'source_sha256': {'foreign': 'b'*64}})
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_branch_parent_path_or_hash_mismatch_rejected(self):
        row = self.row('fine_2')
        row['parent_checkpoint'] = self.rel(self.parent)
        row['parent_sha256'] = self.tx.sha(self.parent)
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_forged_parent_fingerprint_rejected(self):
        self.meta['parent_fingerprint'] = '0'*64
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_nonpositive_nan_or_unbounded_substep_time_rejected(self):
        original = self.row()['dt_s']
        for dt in (0., -1., float('nan'), float('inf'), .5, True):
            with self.subTest(dt=dt):
                self.row()['dt_s'] = dt
                with self.assertRaises((ValueError, AssertionError)):
                    self.validate()
        self.row()['dt_s'] = original

    def test_skipped_physical_step_rejected(self):
        self.row()['step'] += 1
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_changed_old_machine_prefix_rejected_even_with_same_count(self):
        state = self.tx.states[self.row()['sha256']]
        state['machine_history'][0]['dt_s'] = .099
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_changed_old_energy_prefix_rejected_even_with_same_count(self):
        self.tx.states[self.row()['sha256']]['interval_energy']['records'][0] = {'energy': 'replacement'}
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_multiple_appended_records_in_one_leaf_rejected(self):
        state = self.tx.states[self.row()['sha256']]
        state['machine_history'].append(copy.deepcopy(state['machine_history'][-1]))
        state['interval_energy']['records'].append({'extra': True})
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_macro_alias_must_copy_final_leaf_exactly(self):
        alias = self.folder/'coarse/NOT_RESTARTABLE_candidate.npz'
        self.write(alias, b'fabricated-macro-endpoint')
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_macro_observer_must_equal_last_leaf(self):
        self.write(self.folder/'fine_1/observer.json', {'time_s': 99.})
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_final_step_must_count_all_fine_substeps(self):
        self.meta['screened_end_step'] = 4
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_failed_original_pair_gate_rejected(self):
        for part, key in (('global_screen', 'pass'), ('local_screen', 'passed')):
            with self.subTest(part=part):
                self.meta['pair_screen'][part][key] = False
                self.write(self.folder/'pair_screen.json', self.meta['pair_screen'])
                with self.assertRaises((ValueError, AssertionError)):
                    self.validate()
                self.meta['pair_screen'][part][key] = True
        self.meta['pair_screen']['crossing_ids_match'] = False
        self.write(self.folder/'pair_screen.json', self.meta['pair_screen'])
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_event_resolution_bound_to_qualified_leaf(self):
        evidence=self.add_event()
        self.validate()
        self.write(evidence['resolution_path'], {'classification': 'ALTERED_EVENT_RECEIPT_FIXTURE'})
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_unbound_event_resolution_file_rejected(self):
        resolution = (self.root/self.row()['checkpoint']).parent/'event_resolution.json'
        self.write(resolution, {'fixture': True}, artifact=False)
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def test_event_evidence_outside_role_rejected(self):
        self.row()['observer_path'] = self.row('coarse')['observer_path']
        with self.assertRaises((ValueError, AssertionError)):
            self.validate()

    def add_event(self):
        row=self.row();next_row=self.row(index=1)
        high_path=self.root/row['checkpoint'];leaf=high_path.parent
        high_state=self.tx.states[row['sha256']]
        self.tx.particles[row['sha256']]=self.particles(1.+1e-8)
        self.tx.particles[next_row['sha256']]=self.particles(1.+1e-6)
        tolerance=1e-7;low_dt=row['dt_s']-.5e-7
        low_state=copy.deepcopy(high_state)
        low_state['time_s']=row['parent_time_s']+low_dt
        low_state['machine_history'][-1].update(t1_s=low_state['time_s'],dt_s=low_dt)
        low_path=leaf.parent/'event_low/NOT_RESTARTABLE_candidate.npz'
        self.checkpoint(low_path,low_state)
        self.tx.particles[self.tx.sha(low_path)]=self.particles(1.-1e-8)
        qualification=dict(status=receipt.QUALIFIED,accepted=False,restart_allowed=False,delayed_activation_pairs=[[10,20]])
        self.write(low_path.parent/'candidate_qualification.json',qualification)
        width=row['dt_s']-low_dt
        high_x=receipt._positions(self.tx.particles[row['sha256']])
        low_x=receipt._positions(self.tx.particles[self.tx.sha(low_path)])
        parent_pa,parent_state=self.tx.load(self.root/row['parent_checkpoint'])
        delta=float(np.max(np.linalg.norm(high_x-low_x,axis=1)))
        location_path=leaf/'event_location.json'
        location=dict(status='BRACKETED_ACTIVATION',accepted=False,restart_allowed=False,diagnostic_only=True,
            interpolated_states=0,assumed_zero_endpoint_used_as_result=False,
            full_dt_s=.1,time_tolerance_s=tolerance,position_difference_limit_m=5e-8,
            bracket_width_s=width,maximum_position_difference_m=delta,
            lo=dict(blocked=True,dt_s=low_dt,checkpoint=self.rel(low_path),sha256=self.tx.sha(low_path)),
            hi=dict(blocked=False,dt_s=row['dt_s'],checkpoint=row['checkpoint'],sha256=row['sha256']))
        self.write(location_path,location)
        endpoint=1.-tolerance/row['dt_s']
        pose=dict(origin=np.zeros(3),rotation=np.eye(3))
        certificate=audit_blocked_prefix(receipt._positions(parent_pa),high_x,np.array([0]),np.array([1]),
            pose,pose,high_state['parameters'],parent_pa.particle_id,endpoint)
        self.assertEqual(certificate['status'],'CERTIFIED_BLOCKED_PREFIX')
        prefix_path=leaf/'incoming_blocked_prefix.json';self.write(prefix_path,certificate)
        bracket=dict(direction='BLOCKED_TO_VISIBLE',particle_id_i=10,particle_id_j=20,
            start_fraction=1.-.5*tolerance/row['dt_s'],end_fraction=1.)
        routing=dict(status='VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE',unresolved_samples=[],transition_brackets=[bracket])
        self.write(leaf/'union_visibility_events.json',routing)
        self.write(leaf/'visible_path_audit.json',{'status':'CERTIFIED_VISIBLE_SUPPORT_PATH'})
        next_graph=(self.root/next_row['checkpoint']).parent/'visible_graph.npz'
        np.savez(next_graph,particle_id=np.array([10,20]),i=np.array([0]),j=np.array([1]),
            position0=high_x,position1=receipt._positions(self.tx.particles[next_row['sha256']]))
        self.meta['artifacts_sha256'][self.rel(next_graph)]=self.tx.sha(next_graph)
        resolution=dict(status='PASS_BOUNDED_PHYSICALLY_RESOLVED_ACTIVATION',protocol=self.meta['event_protocol'],
            incoming_pid_pairs=[[10,20]],location=self.rel(location_path),location_sha256=self.tx.sha(location_path),
            tolerance_s=tolerance,low_checkpoint=self.rel(low_path),high_checkpoint=row['checkpoint'],
            bracket_width_s=width,real_solution_bracket_width_s=width,maximum_position_difference_m=delta,
            first_post_event_checkpoint=next_row['checkpoint'],incoming_restored=True,original_active_path_certified=True,
            accepted=False,restart_allowed=False,interpolation_accepted=False,
            observed_high_path_activation_delay_s=row['dt_s']*(1.-bracket['start_fraction']),
            activation_delay_bound_s=tolerance,blocked_prefix_path=self.rel(prefix_path),
            blocked_prefix_sha256=self.tx.sha(prefix_path),blocked_prefix_end_fraction=endpoint)
        resolution_path=leaf/'event_resolution.json';self.write(resolution_path,resolution)
        qualification.update(event_resolution_path=self.rel(resolution_path),event_resolution_sha256=self.tx.sha(resolution_path))
        self.write(self.root/row['qualification_path'],qualification)
        self.meta['event_records']['fine_1']=[copy.deepcopy(resolution)]
        return dict(resolution_path=resolution_path,resolution=resolution,location_path=location_path,
            location=location,prefix_path=prefix_path,certificate=certificate,qualification_path=self.root/row['qualification_path'],
            graph_path=next_graph,low_path=low_path)

    def reseal_event(self,evidence):
        self.write(evidence['location_path'],evidence['location'])
        self.write(evidence['prefix_path'],evidence['certificate'])
        resolution=evidence['resolution']
        resolution['location_sha256']=self.tx.sha(evidence['location_path'])
        resolution['blocked_prefix_sha256']=self.tx.sha(evidence['prefix_path'])
        self.write(evidence['resolution_path'],resolution)
        q=self.tx.read(evidence['qualification_path']);q['event_resolution_sha256']=self.tx.sha(evidence['resolution_path'])
        self.write(evidence['qualification_path'],q)
        self.meta['event_records']['fine_1']=[copy.deepcopy(resolution)]

    def test_changed_loaded_leaf_implementation_source_rejected(self):
        self.tx.states[self.row()['sha256']]['implementation_source_sha256']={'foreign':'0'*64}
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_delayed_incoming_leaf_cannot_omit_resolution(self):
        q=self.root/self.row()['qualification_path']
        content=self.tx.read(q);content['delayed_activation_pairs']=[[10,20]];self.write(q,content)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_manifest_event_inventory_must_exactly_match_selected_leaf(self):
        self.add_event();self.meta['event_records']['fine_1']=[]
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_resigned_resolution_cannot_change_selected_high_or_next_leaf(self):
        e=self.add_event();e['resolution']['high_checkpoint']=self.row(index=1)['checkpoint'];self.reseal_event(e)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_resigned_resolution_cannot_claim_interpolation(self):
        e=self.add_event();e['resolution']['interpolation_accepted']=True;self.reseal_event(e)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_resigned_location_must_have_blocked_low_and_visible_high(self):
        e=self.add_event();e['location']['lo']['blocked']=False;self.reseal_event(e)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_actual_endpoint_position_difference_cannot_exceed_bound(self):
        e=self.add_event();self.tx.particles[self.tx.sha(e['low_path'])]=self.particles(1.-1e-5)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_actual_low_visibility_must_be_blocked(self):
        e=self.add_event();self.tx.particles[self.tx.sha(e['low_path'])]=self.particles(1.+.5e-8)
        difference=float(np.max(np.linalg.norm(receipt._positions(self.tx.particles[self.row()['sha256']])-
            receipt._positions(self.tx.particles[self.tx.sha(e['low_path'])]),axis=1)))
        e['location']['maximum_position_difference_m']=difference;e['resolution']['maximum_position_difference_m']=difference
        self.reseal_event(e)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_resigned_prefix_must_bind_actual_high_chord(self):
        e=self.add_event();e['certificate']['input_sha256']['world_position1']='0'*64;self.reseal_event(e)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_resigned_prefix_must_cover_declared_delay(self):
        e=self.add_event();e['certificate']['prefix_end_fraction']=.5
        e['resolution']['blocked_prefix_end_fraction']=.5;self.reseal_event(e)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_resigned_prefix_cannot_have_missing_certificate_coverage(self):
        e=self.add_event();e['certificate']['certificates']=[];self.reseal_event(e)
        with self.assertRaises((ValueError,AssertionError)):self.validate()

    def test_resigned_post_event_graph_must_restore_incoming_pair(self):
        e=self.add_event()
        with np.load(e['graph_path'],allow_pickle=False) as graph:
            data={key:graph[key] for key in graph.files}
        data['i']=np.array([],dtype=int);data['j']=np.array([],dtype=int)
        np.savez(e['graph_path'],**data)
        self.meta['artifacts_sha256'][self.rel(e['graph_path'])]=self.tx.sha(e['graph_path'])
        with self.assertRaises((ValueError,AssertionError)):self.validate()


if __name__ == '__main__':
    unittest.main(verbosity=2)
