"""Real Windows sharing-lock tests; only temporary manufactured files are used."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from afsd_locked_history_20261010 import LockedHistory, HistoryLockError


MUTATE = r'''
from pathlib import Path
import json, os, sys
op, name, target = sys.argv[1:]
p = Path(name)
try:
    if op == 'write':
        p.write_bytes(b'CHANGED!')
    elif op == 'delete':
        p.unlink()
    elif op == 'replace':
        os.replace(target, p)
    elif op == 'rename':
        p.rename(target)
    elif op == 'create':
        p.mkdir()
        (p / 'new.bin').write_bytes(b'new output')
    else:
        raise RuntimeError(op)
except OSError as error:
    print(json.dumps({'blocked': True, 'winerror': error.winerror,
                      'errno': error.errno, 'error': str(error)}))
else:
    print(json.dumps({'blocked': False}))
'''


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


@unittest.skipUnless(os.name == 'nt', 'Windows locking guarantees require Windows')
class LockedHistoryTests(unittest.TestCase):
    def setUp(self):
        # Use the actual project volume. The app sandbox virtualizes parts of
        # the user profile and does not allow opening every ancestor of %TEMP%.
        self.temp = tempfile.TemporaryDirectory(prefix='afsd_locked_history_',
                                                dir=ROOT / 'outputs')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.parent = self.root / 'evidence'
        self.parent.mkdir()
        self.file = self.parent / 'file.bin'
        self.payload = b'ORIGINAL'
        self.file.write_bytes(self.payload)
        self.digest = sha(self.payload)

    def pool(self):
        pool = LockedHistory(self.root)
        self.addCleanup(pool.close)
        return pool

    def mutate(self, op, path=None, target=''):
        child = subprocess.run([sys.executable, '-B', '-c', MUTATE, op,
                                str(path or self.file), str(target)],
                               capture_output=True, text=True, timeout=15, check=True)
        return json.loads(child.stdout)

    def test_cache_reuse_requires_live_identity_without_rereading_bytes(self):
        pool = self.pool()
        pool.seal({'evidence/file.bin': self.digest})
        first = pool.stats()
        pool.verify({str(self.file): self.digest.upper()})
        self.assertEqual(pool.sha(self.file), self.digest)
        later = pool.stats()
        self.assertEqual(later['hash_reads'], 1)
        self.assertEqual(later['hash_bytes'], len(self.payload))
        self.assertEqual(later['unique_files'], 1)
        self.assertEqual(later['cache_hits'], 2)
        self.assertGreater(later['identity_checks'], first['identity_checks'])

    def test_same_size_write_blocked_in_separate_process_and_close_releases(self):
        pool = self.pool()
        pool.seal({self.file: self.digest})
        result = self.mutate('write')
        self.assertTrue(result['blocked'])
        self.assertTrue(result['winerror'] in (5, 32, 33) or result['errno'] == 13)
        self.assertEqual(self.file.read_bytes(), self.payload)
        pool.verify({self.file: self.digest})
        pool.close()
        self.assertFalse(self.mutate('write')['blocked'])
        self.assertEqual(self.file.read_bytes(), b'CHANGED!')

    def test_delete_blocked_in_separate_process(self):
        pool = self.pool()
        pool.seal({self.file: self.digest})
        self.assertTrue(self.mutate('delete')['blocked'])
        self.assertTrue(self.file.exists())

    def test_atomic_replace_blocked_in_separate_process(self):
        replacement = self.parent / 'replacement.bin'
        replacement.write_bytes(b'CHANGED!')
        pool = self.pool()
        pool.seal({self.file: self.digest})
        self.assertTrue(self.mutate('replace', target=replacement)['blocked'])
        self.assertEqual(self.file.read_bytes(), self.payload)
        self.assertTrue(replacement.exists())

    def test_ancestor_rename_blocked_and_new_outputs_allowed(self):
        pool = self.pool()
        pool.seal({self.file: self.digest})
        self.assertTrue(self.mutate('rename', self.parent,
                                    self.root / 'moved')['blocked'])
        self.assertTrue(self.mutate('rename', self.root,
                                    str(self.root) + '_moved')['blocked'])
        self.assertFalse(self.mutate('create', self.parent / 'new_stage')['blocked'])
        self.assertEqual((self.parent / 'new_stage/new.bin').read_bytes(), b'new output')

    def test_live_write_handle_prevents_sealing(self):
        pool = self.pool()
        with self.file.open('r+b'):
            with self.assertRaises(OSError):
                pool.seal({self.file: self.digest})
        self.assertTrue(pool.stats()['failed'])
        self.assertEqual(pool.stats()['locked_files'], 0)

    def test_hash_mismatch_invalidates_session_and_releases_all_files(self):
        second = self.parent / 'other.bin'
        second.write_bytes(b'other')
        pool = self.pool()
        pool.seal({self.file: self.digest})
        with self.assertRaisesRegex(HistoryLockError, 'hash mismatch'):
            pool.verify({second: '0' * 64})
        self.assertTrue(pool.stats()['failed'])
        self.assertEqual(pool.stats()['locked_files'], 0)
        self.assertFalse(self.mutate('write')['blocked'])
        with self.assertRaisesRegex(HistoryLockError, 'closed or invalidated'):
            pool.sha(self.file)

    def test_conflicting_expectation_rejected_before_cached_authorization(self):
        pool = self.pool()
        pool.seal({self.file: self.digest})
        with self.assertRaisesRegex(HistoryLockError, 'conflicting expected'):
            pool.verify({'evidence/file.bin': 'f' * 64})
        self.assertTrue(pool.stats()['closed'])

    def test_case_alias_uses_single_file_and_conflict_is_rejected(self):
        pool = self.pool()
        pool.seal({str(self.file): self.digest})
        pool.verify({str(self.file).upper(): self.digest})
        self.assertEqual(pool.stats()['hash_reads'], 1)
        with self.assertRaisesRegex(HistoryLockError, 'conflicting expected'):
            pool.verify({str(self.file).upper(): 'f' * 64})

    def test_sha_then_wrong_expected_digest_rejected(self):
        pool = self.pool()
        self.assertEqual(pool.sha(self.file), self.digest)
        with self.assertRaisesRegex(HistoryLockError, 'hash mismatch'):
            pool.verify({self.file: '0' * 64})

    def test_final_audit_rereads_all_bytes_and_retains_locks(self):
        large = self.parent / 'large.bin'
        payload = bytes(range(256)) * 8193
        large.write_bytes(payload)
        pool = self.pool()
        pool.seal({self.file: self.digest, large: sha(payload)})
        pool.final_audit()
        stats = pool.stats()
        self.assertEqual(stats['hash_reads'], 4)
        self.assertEqual(stats['hash_bytes'], 2 * (len(self.payload) + len(payload)))
        self.assertEqual(stats['final_audits'], 1)
        self.assertEqual(stats['locked_files'], 2)
        self.assertTrue(self.mutate('write')['blocked'])

    def test_read_bytes_uses_held_handle_and_matches_digest(self):
        pool = self.pool()
        pool.seal({self.file: self.digest})
        self.assertEqual(pool.read_bytes(self.file), self.payload)
        self.assertEqual(pool.stats()['hash_reads'], 2)

    def test_reopened_different_file_identity_invalidates_cache(self):
        # Real filesystem replacement is denied by earlier tests. Inject only
        # name resolution here, returning a real handle to another same-sized
        # file, to verify the independent identity guard also fails closed.
        other = self.parent / 'other.bin'
        other.write_bytes(self.payload)
        pool = self.pool()
        pool.seal({self.file: self.digest})
        native_open = pool._win.open
        with mock.patch.object(pool._win, 'open',
                               side_effect=lambda path, directory=False:
                               native_open(other if path == self.file else path, directory)):
            with self.assertRaisesRegex(HistoryLockError, 'different file'):
                pool.sha(self.file)
        self.assertTrue(pool.stats()['failed'])

    def test_dead_held_handle_invalidates_cache(self):
        pool = self.pool()
        pool.seal({self.file: self.digest})
        row = next(iter(pool._files.values()))
        pool._win.close(row.handle)
        with self.assertRaises(OSError):
            pool.sha(self.file)
        self.assertTrue(pool.stats()['failed'])

    def test_final_audit_detects_digest_corruption(self):
        pool = self.pool()
        pool.seal({self.file: self.digest})
        next(iter(pool._files.values())).digest = '0' * 64
        with self.assertRaisesRegex(HistoryLockError, 'final byte audit mismatch'):
            pool.final_audit()

    def test_noncanonical_paths_and_escape_are_rejected(self):
        paths = ['../outside.bin', 'evidence/../evidence/file.bin',
                 './evidence/file.bin', 'evidence/file.bin:stream',
                 'evidence/file.bin.', 'evidence/file.bin ',
                 str(self.root.parent / 'outside.bin'), 'E:relative.bin',
                 '\\\\?\\' + str(self.file), str(self.root)]
        for name in paths:
            with self.subTest(name=name):
                pool = self.pool()
                with self.assertRaises(HistoryLockError):
                    pool.sha(name)
                self.assertTrue(pool.stats()['closed'])

    def test_invalid_hash_and_missing_file_fail_closed(self):
        for digest in ('x' * 64, '', 23, 'a' * 63):
            with self.subTest(digest=digest):
                pool = self.pool()
                with self.assertRaisesRegex(HistoryLockError, 'valid SHA256'):
                    pool.verify({self.file: digest})
        pool = self.pool()
        with self.assertRaises(OSError):
            pool.sha('missing.bin')
        self.assertTrue(pool.stats()['closed'])

    def test_directory_is_not_a_file(self):
        pool = self.pool()
        with self.assertRaises((OSError, HistoryLockError)):
            pool.sha(self.parent)
        self.assertTrue(pool.stats()['failed'])

    def test_context_exit_releases_handles_and_old_cache_is_unusable(self):
        with LockedHistory(self.root) as pool:
            pool.seal({self.file: self.digest})
            self.assertTrue(self.mutate('write')['blocked'])
        self.assertFalse(self.mutate('write')['blocked'])
        with self.assertRaisesRegex(HistoryLockError, 'closed or invalidated'):
            pool.verify({self.file: self.digest})

    def test_process_exit_releases_locks_without_python_cleanup(self):
        code = r'''
import os, sys
sys.path.insert(0, sys.argv[1])
from afsd_locked_history_20261010 import LockedHistory
pool = LockedHistory(sys.argv[2])
pool.seal({sys.argv[3]: sys.argv[4]})
print('READY', flush=True)
sys.stdin.readline()
os._exit(0)
'''
        child = subprocess.Popen([sys.executable, '-B', '-c', code,
                                  str(ROOT / 'scripts'), str(self.root),
                                  str(self.file), self.digest],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, text=True,
                                 creationflags=subprocess.CREATE_NO_WINDOW)
        try:
            self.assertEqual(child.stdout.readline().strip(), 'READY')
            self.assertTrue(self.mutate('write')['blocked'])
            _, error = child.communicate('\n', timeout=15)
            self.assertEqual(child.returncode, 0, error)
            self.assertFalse(self.mutate('write')['blocked'])
        finally:
            if child.poll() is None:
                child.kill()
                child.communicate(timeout=15)

    def test_existing_hard_link_cannot_bypass_write_lock(self):
        alias = self.root / 'hardlink.bin'
        os.link(self.file, alias)
        pool = self.pool()
        pool.seal({self.file: self.digest})
        self.assertTrue(self.mutate('write', alias)['blocked'])
        self.assertEqual(alias.read_bytes(), self.payload)

    def test_directory_symlink_is_rejected_when_creation_is_available(self):
        alias = self.root / 'symlink'
        try:
            alias.symlink_to(self.parent, target_is_directory=True)
        except OSError as error:
            if error.winerror == 1314:
                self.skipTest('Windows symlink creation privilege unavailable')
            raise
        pool = self.pool()
        with self.assertRaisesRegex(HistoryLockError, 'reparse points'):
            pool.sha(alias / self.file.name)

    def test_ten_thousand_native_file_handles_and_no_hash_reads_on_reverify(self):
        # Manufactured one-byte fixtures test the actual Windows handle count;
        # no project history or physical simulation data are opened here.
        many = self.root / 'many'
        many.mkdir()
        mapping = {}
        digest = sha(b'x')
        for index in range(10003):
            file = many / ('%05d.bin' % index)
            file.write_bytes(b'x')
            mapping[file] = digest
        pool = self.pool()
        pool.seal(mapping)
        self.assertEqual(pool.stats()['locked_files'], 10003)
        self.assertEqual(pool.stats()['hash_reads'], 10003)
        pool.verify(mapping)
        self.assertEqual(pool.stats()['hash_reads'], 10003)
        self.assertEqual(pool.stats()['cache_hits'], 10003)
        pool.final_audit()
        self.assertEqual(pool.stats()['hash_reads'], 20006)
        self.assertEqual(pool.stats()['hash_bytes'], 20006)


if __name__ == '__main__':
    unittest.main(verbosity=2)
