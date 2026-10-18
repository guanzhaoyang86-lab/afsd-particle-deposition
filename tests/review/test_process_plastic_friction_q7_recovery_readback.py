"""SCHEMA_FIXTURE_NOT_PHYSICAL: Q7R metadata gate, mocked parent physics audit.

These tests do not create or validate a material/contact solution. The frozen
parent saved-data audit is deliberately mocked to isolate rejection ordering.
"""
import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import inspect_process_plastic_friction_recovery as reader


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


class RecoveryMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name)
        self.recovery = reader.read_json(ROOT/reader.RECOVERY_CONFIG)
        self.cfg = reader.read_json(ROOT/reader.BASELINE)
        self.cfg['plastic_parent']['numerical']['maximum_iterations'] = 14
        self.cfg['recovery_variant'] = copy.deepcopy(self.recovery)
        self.sources = {name:reader.digest(ROOT/name) for name in reader.BASELINE_SOURCES}
        self.sidecar = dict(variant=self.recovery['variant'], production_authorized=False,
            configuration=copy.deepcopy(self.recovery), baseline_source_sha256=self.sources.copy(),
            source_sha256={name:reader.digest(ROOT/name) for name in reader.SUPPLEMENTAL_SOURCES})
        self.report = dict(increments_per_phase=4, source_sha256=self.sources.copy(), artifact_sha256={})
        self.write()
        self.delegate = patch.object(reader.parent, 'inspect', return_value={
            'status':'PASS_COMPLETE_READBACK', 'evidence_kind':'SCHEMA_FIXTURE_NOT_PHYSICAL'})
        self.mock_parent = self.delegate.start()
        self.addCleanup(self.delegate.stop)

    def write(self):
        save(self.path/'resolved_config.json', self.cfg)
        save(self.path/'recovery_provenance.json', self.sidecar)
        self.report['artifact_sha256']['recovery_provenance.json'] = reader.digest(self.path/'recovery_provenance.json')
        save(self.path/'report.json', self.report)

    def rejected(self, pattern):
        with self.assertRaisesRegex(ValueError, pattern):
            reader.inspect(self.path)
        self.mock_parent.assert_not_called()

    def test_complete_metadata_delegates_without_solver_reference(self):
        result = reader.inspect(self.path)
        self.assertEqual(result['status'], 'PASS_COMPLETE_READBACK')
        self.assertEqual(result['recovery_audit']['maximum_iterations'], 14)
        self.assertTrue(result['recovery_audit']['physics_tolerances_boundaries_unchanged'])
        self.mock_parent.assert_called_once_with(self.path)

    def test_eight_increment_grid_is_only_other_allowed_case(self):
        self.cfg['increments_per_phase'] = 8
        self.report['increments_per_phase'] = 8
        self.write()
        self.assertEqual(reader.inspect(self.path)['status'], 'PASS_COMPLETE_READBACK')

    def test_physical_change_rejected_even_with_rehashed_metadata(self):
        self.cfg['parameters']['young']['value'] *= .9
        self.write()
        self.rejected('configuration changed outside')

    def test_nested_contact_tolerance_change_rejected(self):
        self.cfg['friction_parent']['numerical']['force_tolerance_N'] *= 10
        self.write()
        self.rejected('configuration changed outside')

    def test_unapproved_newton_budget_rejected(self):
        self.cfg['plastic_parent']['numerical']['maximum_iterations'] = 15
        self.write()
        self.rejected('configuration changed outside')

    def test_recovery_metadata_missing_from_resolved_config_rejected(self):
        del self.cfg['recovery_variant']
        self.write()
        self.rejected('configuration changed outside')

    def test_supplemental_entry_omission_rejected(self):
        self.sidecar['source_sha256'].pop(next(iter(self.sidecar['source_sha256'])))
        self.write()
        self.rejected('supplemental source inventory')

    def test_baseline_digest_tamper_rejected(self):
        self.sidecar['baseline_source_sha256'][next(iter(self.sources))] = '0'*64
        self.write()
        self.rejected('baseline source inventory')

    def test_recovery_provenance_uncovered_by_report_rejected(self):
        self.report['artifact_sha256'] = {}
        save(self.path/'report.json', self.report)
        self.rejected('missing from artifact inventory')

    def test_post_manifest_sidecar_tamper_rejected(self):
        self.sidecar['note'] = 'added after report digest was recorded'
        save(self.path/'recovery_provenance.json', self.sidecar)
        self.rejected('missing from artifact inventory or changed')

    def test_recovery_variant_metadata_mismatch_rejected(self):
        self.sidecar['variant'] = 'UNDECLARED'
        self.write()
        self.rejected('provenance/configuration mismatch')

    def test_partial_cli_exit_two_complete_exit_zero_and_no_overwrite(self):
        self.mock_parent.return_value = dict(status='PASS_PARTIAL_READBACK',
            evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
        partial = self.path/'partial.json'
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(reader.main(['--run',str(self.path),'--output',str(partial)]),2)
        self.assertEqual(reader.read_json(partial)['status'],'PASS_PARTIAL_READBACK')
        self.mock_parent.return_value = dict(status='PASS_COMPLETE_READBACK',
            evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(reader.main(['--run',str(self.path),'--output',str(self.path/'complete.json')]),0)
        with self.assertRaises(FileExistsError):
            reader.main(['--run',str(self.path),'--output',str(partial)])


if __name__ == '__main__':
    unittest.main()
