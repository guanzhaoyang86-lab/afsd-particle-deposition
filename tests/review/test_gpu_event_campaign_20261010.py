"""Event controller integration/registration checks; no GPU campaign launches."""
import argparse
import ast
import copy
from datetime import datetime, timedelta, timezone
import hashlib
import inspect
from pathlib import Path
import sys
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_event_campaign_20261010 as driver
import afsd_event_worker_20261010 as worker
import check_gpu_visibility_event_20261010 as checker


class EventCampaignTests(unittest.TestCase):
    def fixture(self, parent=None):
        old = {'old-physical-source': 'a'*64}
        new = dict(old, **{'new-event-physical-source': 'b'*64})
        parent = parent or driver.REGISTERED_OLD_REPLAY_PARENT
        meta = dict(source_sha256=old, search_tensor_protocol=driver.search_bridge.PROTOCOL,
            screened_continuation_checkpoint=driver.rel(parent/'fine_2/NOT_RESTARTABLE_candidate.npz'),
            screened_continuation_sha256='c'*64, maximum_candidate_wall_s=40., dt_s=.000474609375,
            screened_end_step=480)
        template = dict(source_sha256=new, search_source_sha256=old, event_components_sha256='d'*64,
            event_tolerance_scale=.5)
        return parent, meta, template

    def test_original_live_supervision_and_outer_receipt_checks_are_preserved(self):
        report = driver.static_audit()
        self.assertTrue(report['old_outer_receipt_checks_unchanged'])
        self.assertTrue(report['live_supervision_unchanged'])
        self.assertFalse(report['runtime_global_rebinding'])

    def test_worker_only_rebinds_explicit_new_provider_and_pair_entrypoint(self):
        source = (ROOT/'scripts/afsd_incremental_worker_resume_20261010.py').read_text(encoding='utf-8')
        source = source.replace('Isolated unchanged physical pair; live parent grants replace history rereads.',
                                'Isolated explicit event pair; original live-history authorization and limits.')
        source = source.replace('import run_gpu_search_tensor_20261009 as physical',
                                'import run_gpu_visibility_event_20261010 as physical')
        source = source.replace('physical._run_pair_original(output,protocol,report,began+525.)',
                                'physical.run_pair(output,protocol,report,began+525.)')
        actual = (ROOT/'scripts/afsd_event_worker_20261010.py').read_text(encoding='utf-8')
        self.assertEqual(ast.dump(ast.parse(source)), ast.dump(ast.parse(actual)))

    def test_old42_to_new_registration_keeps_both_sources_and_exact_parent(self):
        parent, meta, template = self.fixture()
        protocol = driver.stage_protocol(template, meta, parent, meta['dt_s'], {}, {}, {})
        self.assertEqual(protocol['source_sha256'], template['source_sha256'])
        self.assertEqual(protocol['parent_implementation_source_sha256'], meta['source_sha256'])
        transition = protocol['implementation_transition']
        self.assertEqual(transition['from_source_sha256'], meta['source_sha256'])
        self.assertEqual(transition['to_source_sha256'], template['source_sha256'])
        self.assertEqual(transition['parent_sha256'], meta['screened_continuation_sha256'])
        self.assertEqual(transition['components_sha256'], template['event_components_sha256'])
        self.assertEqual(protocol['event_tolerance_scale'], .5)
        self.assertEqual(protocol['event_protocol'], driver.PROTOCOL)

    def test_unregistered_old_parent_cannot_enter_new_implementation(self):
        parent, meta, template = self.fixture(driver.BASE/'campaign_10h_20261010_v1/attempt_0044')
        with self.assertRaises(ValueError):
            driver.source_transition(template, meta, parent)

    def test_foreign_source_or_missing_new_protocol_rejected(self):
        parent, meta, template = self.fixture()
        meta['source_sha256'] = {'foreign-source': 'e'*64}
        with self.assertRaises(ValueError):
            driver.source_transition(template, meta, parent)
        meta['source_sha256'] = template['source_sha256']
        with self.assertRaises(ValueError):
            driver.source_transition(template, meta, parent)
        meta['event_protocol'] = driver.PROTOCOL
        self.assertIsNone(driver.source_transition(template, meta, parent))

    def test_source_closure_never_overwrites_historical_digest(self):
        parent, meta, template = self.fixture()
        with self.assertRaises(ValueError):
            driver.stage_protocol(template, meta, parent, meta['dt_s'], {'old-physical-source': '0'*64}, {}, {})

    def test_three_registered_generations_and_protocol_identity(self):
        old = {'base': 'a'*64}; search = dict(old, search='b'*64); new = dict(search, event='c'*64)
        self.assertEqual(driver.source_generation({'source_sha256':old}, old, search, new), 0)
        self.assertEqual(driver.source_generation({'source_sha256':search,
            'search_tensor_protocol':driver.search_bridge.PROTOCOL}, old, search, new), 1)
        self.assertEqual(driver.source_generation({'source_sha256':new, 'event_protocol':driver.PROTOCOL,
            'event_substeps':{}}, old, search, new), 2)
        for meta in ({'source_sha256':search}, {'source_sha256':new}, {'source_sha256':{'other':'e'*64}}):
            with self.subTest(meta=meta), self.assertRaises(ValueError):
                driver.source_generation(meta, old, search, new)

    def test_forecast_limit_450_seconds_unchanged(self):
        parent, meta, template = self.fixture()
        meta['maximum_candidate_wall_s'] = (450.-30.)/4.5
        self.assertEqual(driver.stage_protocol(template, meta, parent, meta['dt_s'], {}, {}, {})['forecast_s'],450.)
        meta['maximum_candidate_wall_s'] += .001
        with self.assertRaises(ValueError):
            driver.stage_protocol(template, meta, parent, meta['dt_s'], {}, {}, {})

    def test_original_absolute_deadline_cannot_be_extended_by_restart(self):
        now = driver.ABSOLUTE_STOP-timedelta(hours=2)
        self.assertEqual(driver.authorized_wall_seconds(36000.,now),7200.)
        self.assertEqual(driver.authorized_wall_seconds(120.,now),120.)
        for seconds in (0.,-1.,36000.001,float('inf'),float('nan')):
            with self.subTest(seconds=seconds), self.assertRaises(ValueError):
                driver.authorized_wall_seconds(seconds,now)
        with self.assertRaises(ValueError):
            driver.authorized_wall_seconds(1.,driver.ABSOLUTE_STOP)
        with self.assertRaises(ValueError):
            driver.authorized_wall_seconds(1.,datetime(2026,10,10))

    def test_cli_single_stage_and_registered_tolerance_scales(self):
        for scale in ('1', '.5', '.25'):
            argv = ['event-driver','--output','unused','--parent','unused-parent',
                    '--max-stages','1','--event-tolerance-scale',scale]
            with patch.object(sys,'argv',argv), patch.object(driver,'run',return_value=0) as run:
                self.assertEqual(driver.main(),0)
            self.assertEqual(run.call_args.args[0].max_stages,1)
            self.assertEqual(run.call_args.args[0].event_tolerance_scale,float(scale))

    def test_cli_rejects_unregistered_tolerance_scale(self):
        argv = ['event-driver','--output','unused','--parent','unused-parent','--event-tolerance-scale','2']
        with patch.object(sys,'argv',argv), patch('sys.stderr',new=Mock()), self.assertRaises(SystemExit) as cm:
            driver.main()
        self.assertEqual(cm.exception.code,2)

    def test_transaction_root_and_original_outer_receipt_method_are_inherited(self):
        transaction = driver.LockedPreflight(Mock())
        self.assertEqual(transaction.root,driver.ROOT)
        self.assertIs(driver.LockedPreflight._gpu_receipt_uncached,
                      driver.cache.CachedPreflight._gpu_receipt_uncached)

    def test_observation_visits_every_selected_fine_leaf_and_no_coarse_leaf(self):
        rows = {role:[{'checkpoint':role+str(i), 'observer_path':role+str(i)+'.json'} for i in (1,2)]
                for role in ('coarse','fine_1','fine_2')}
        meta = {'event_substeps':{'schema':driver.event_receipt.SCHEMA,**rows}}
        frames = ['initial']; seen = []
        with patch.object(driver,'read',side_effect=lambda p:{'from':p.name}), \
             patch.object(driver,'load',side_effect=lambda p:(p.name, {'parameters':{}})), \
             patch.object(driver.deposition,'frame_from_saved',side_effect=lambda pa,*args:seen.append(pa) or pa), \
             patch.object(driver.deposition,'summarize_progress',side_effect=lambda fs,p:{'frames':list(fs)}):
            report,state = driver.observe(meta,frames)
        self.assertEqual(seen,['fine_11','fine_12','fine_21','fine_22'])
        self.assertEqual(report['frames'],['initial']+seen)

    def test_validated_stage_counts_actual_fine_physical_steps(self):
        parent,meta,template=self.fixture()
        stage = dict(source_sha256=template['source_sha256'],parent_sha256=meta['screened_continuation_sha256'],
            parent_implementation_source_sha256=meta['source_sha256'],event_protocol=driver.PROTOCOL,
            event_substeps={'fine_1':[{},{}],'fine_2':[{}, {},{}]},screened_end_step=485)
        tx=Mock();tx.gpu_receipt.return_value=stage
        with patch.object(driver,'LockedPreflight',return_value=tx):
            self.assertEqual(driver.validated_stage(Mock(),driver.BASE/'test-stage',template['source_sha256'],meta),stage)
            stage['screened_end_step']=482
            with self.assertRaises(ValueError):
                driver.validated_stage(Mock(),driver.BASE/'test-stage',template['source_sha256'],meta)

    def geometry_fixture(self, attempt='incremental_deposition_resume_20261010_v1/attempt_0019'):
        folder=driver.BASE/attempt
        failure=driver.read(folder/'manifest.json')
        parent=driver.read(driver.ROOT/failure['parent_gpu_directory']/'manifest.json')
        policy=driver.adaptive.initial_state(parent['screened_continuation_sha256'],failure['dt_s'])
        policy=driver.adaptive.begin_attempt(policy,failure['dt_s'])
        pool=Mock();pool.sha.side_effect=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
        return folder,failure,parent,policy,pool

    def test_original_pure_budget_failure_remains_same_parent_refinable(self):
        folder,failure,parent,policy,pool=self.geometry_fixture()
        smaller,event=driver.geometry_budget_recovery(pool,folder,failure,parent,policy,parent['source_sha256'])
        self.assertEqual(event['action'],'REFINE')
        self.assertEqual(smaller['parent_sha256'],parent['screened_continuation_sha256'])
        self.assertEqual(smaller['nominal_dt_s'],failure['dt_s']/2)

    def test_actual_old44_activation_is_not_misclassified_as_pure_budget(self):
        folder,failure,parent,policy,pool=self.geometry_fixture('campaign_10h_20261010_v1/attempt_0044')
        self.assertIsNone(driver.geometry_budget_recovery(pool,folder,failure,parent,policy,parent['source_sha256']))

    def test_failed_geometry_must_use_expected_new_source(self):
        folder,failure,parent,policy,pool=self.geometry_fixture()
        with self.assertRaises(ValueError):
            driver.geometry_budget_recovery(pool,folder,failure,parent,policy,{'foreign':'0'*64})

    def test_checker_and_controller_physical_registration_match(self):
        self.assertEqual(set(checker.PHYSICAL_FILES),set(driver.PHYSICAL_FILES))
        self.assertIn('tests/review/test_gpu_event_campaign_20261010.py',checker.TEST_FILES)
        self.assertIn('tests/review/test_gpu_visibility_event_20261010.py',checker.TEST_FILES)

    def preflight_fixture(self):
        legacy={'legacy':'a'*64}; search=dict(legacy,search='b'*64); current=dict(search,event='c'*64)
        valid_dir=driver.BASE/'replay_v1'
        old_dir=driver.REGISTERED_OLD_REPLAY_PARENT
        new_dir=driver.BASE/'unexecuted_event_fixture/attempt_0001'
        def meta(source,folder,digest):
            return dict(source_sha256=source,validation_directory=driver.rel(valid_dir),
                screened_continuation_checkpoint=driver.rel(folder/'fine_2/NOT_RESTARTABLE_candidate.npz'),
                screened_continuation_sha256=digest,maximum_candidate_wall_s=40.,dt_s=.000474609375)
        valid=meta(legacy,valid_dir,'1'*64)
        old=meta(search,old_dir,'2'*64)
        old.update(parent_sha256=valid['screened_continuation_sha256'],parent_gpu_directory=None,
            parent_implementation_source_sha256=legacy,search_tensor_protocol=driver.search_bridge.PROTOCOL)
        new=meta(current,new_dir,'3'*64)
        new.update(parent_sha256=old['screened_continuation_sha256'],parent_gpu_directory=driver.rel(old_dir),
            parent_implementation_source_sha256=search,event_protocol=driver.PROTOCOL,event_substeps={})
        template=dict(source_sha256=current,search_source_sha256=search,event_components_sha256='4'*64)
        new['implementation_transition']=driver.source_transition(template,old,old_dir)
        records={valid_dir.resolve():valid,old_dir.resolve():old,new_dir.resolve():new}
        transaction=Mock()
        transaction.preflight.return_value=dict(source_sha256=legacy,reference_sha256={})
        transaction.gpu_receipt.side_effect=lambda folder,statuses:copy.deepcopy(records[Path(folder).resolve()])
        transaction.sha.return_value='4'*64
        transaction.finish.return_value={'status':'SYNTHETIC_PREFLIGHT'}
        args=argparse.Namespace(mode='continue',parent=driver.rel(new_dir),dt=None,
            validation=driver.rel(valid_dir),reference_stage=None)
        return args,transaction,records,legacy,search,current,new_dir,old_dir

    def execute_preflight_fixture(self, fixture):
        args,tx,records,legacy,search,current,new_dir,old_dir=fixture
        with patch.object(driver,'LockedPreflight',return_value=tx), \
             patch.object(driver.search_bridge,'frozen_sources',return_value=search), \
             patch.object(driver,'event_sources',return_value=current), \
             patch.object(driver.cache,'bind_administration'), \
             patch.object(driver,'adaptation_audit',return_value={}), \
             patch.object(driver.search_bridge,'adaptation_audit',return_value={}):
            return driver._initial_preflight(args,Mock())

    def test_preflight_preserves_old_to_new_chain_and_exact_initial_parent(self):
        fixture=self.preflight_fixture()
        result=self.execute_preflight_fixture(fixture)
        args,tx,records,legacy,search,current,new_dir,old_dir=fixture
        self.assertEqual(result['source_sha256'],current)
        self.assertEqual(result['legacy_source_sha256'],legacy)
        self.assertEqual(result['search_source_sha256'],search)
        self.assertEqual(result['parent_sha256'],records[new_dir]['screened_continuation_sha256'])
        self.assertEqual(result['event_components_sha256'],'4'*64)
        tx.finish.assert_called_once()
        self.assertEqual(tx.references.call_count,2)

    def test_preflight_rejects_changed_parent_source_and_boundary_receipt(self):
        for kind in ('parent_sha','parent_source','transition_sha','transition_from','event_protocol'):
            fixture=self.preflight_fixture()
            args,tx,records,legacy,search,current,new_dir,old_dir=fixture
            new=records[new_dir]
            if kind=='parent_sha':new['parent_sha256']='0'*64
            elif kind=='parent_source':new['parent_implementation_source_sha256']=legacy
            elif kind=='transition_sha':new['implementation_transition']['components_sha256']='0'*64
            elif kind=='transition_from':new['implementation_transition']['from_source_sha256']=legacy
            else:new['event_protocol']='unregistered'
            with self.subTest(kind=kind),self.assertRaises(ValueError):
                self.execute_preflight_fixture(fixture)

    def test_preflight_rejects_reverting_from_event_generation_to_search_generation(self):
        fixture=self.preflight_fixture()
        args,tx,records,legacy,search,current,new_dir,old_dir=fixture
        # Newest entry claims search while its predecessor is already event.
        records[new_dir].update(source_sha256=search,search_tensor_protocol=driver.search_bridge.PROTOCOL,
            parent_implementation_source_sha256=current)
        records[old_dir].update(source_sha256=current,event_protocol=driver.PROTOCOL,event_substeps={})
        with self.assertRaises(ValueError):
            self.execute_preflight_fixture(fixture)

    def test_event_source_registration_requires_every_physical_file_and_no_conflicts(self):
        tx=Mock();tx.sha.return_value='9'*64
        physical={name:'1'*64 for name in driver.PHYSICAL_FILES}
        components=dict(status='PASS_INCREMENTAL_COMPONENTS',exit_code=0,physical_source_sha256=physical,
            source_sha256=dict(physical),artifacts_sha256={})
        tx.read.return_value=components
        result=driver.event_sources(tx,{'old':'a'*64})
        self.assertEqual(result['old'],'a'*64)
        self.assertEqual(result[driver.rel(driver.INCREMENTAL_CHECKS)],'9'*64)
        self.assertEqual(set(physical)-set(result),set())
        components['source_sha256'][driver.PHYSICAL_FILES[0]]='0'*64
        with self.assertRaises(ValueError):driver.event_sources(tx,{'old':'a'*64})


if __name__ == '__main__':
    unittest.main(verbosity=2)
