"""Per-launch deduplication of administrative evidence checks only.

Each new launch hashes every unique dependency from bytes, retains every old
receipt/history assertion, and rereads every dependency at the final seal.
No persistent cache or stat-only authorization exists. Physics runs exclusively
in the unchanged run_afsd_gpu_20261007.py --worker subprocess.
"""
from pathlib import Path
import argparse
import ast
import copy
import hashlib
import inspect
import json
import os
import subprocess
import sys
import textwrap
import time
import numpy as np
import run_afsd_gpu_20261007 as original

ROOT, BASE, CPU_BASE = original.ROOT, original.BASE, original.CPU_BASE
cpu, core, model, gpu = original.cpu, original.core, original.model, original.gpu
require, rel = original.require, original.rel
CONTRACT = original.CONTRACT
BENCH_PASS, REPLAY_PASS, CONTINUE_PASS = original.BENCH_PASS, original.REPLAY_PASS, original.CONTINUE_PASS
COST_PASS, STAGE_PASS, SHOULDER_PASS, PREFIX_PASS = cpu.COST_PASS, cpu.STAGE_PASS, cpu.SHOULDER_PASS, cpu.PREFIX_PASS
ORIGINAL_RUNNER = Path(original.__file__).resolve()
ADMIN_CONTRACT = ROOT/'docs/GPU_FAST_CONTINUATION_20261008.md'
DEFAULT_CHECKS = BASE/'preflight_components_20261008_v1.json'
adaptation_audit = original.adaptation_audit


class CachedPreflight:
    """A disposable read-only transaction; cached loads return detached copies."""
    def __init__(self, root=ROOT):
        self.root = Path(root).resolve()
        self.records = {}
        self.expected = {}
        self.json_cache = {}
        self.load_cache = {}
        self.cpu_receipts = {}
        self.gpu_receipts = {}
        self.requests = 0
        self.hash_reads = 0
        self.hash_bytes = 0
        self.finished = False
        self.started = time.perf_counter()

    def path(self, value):
        path = Path(value)
        path = (self.root/path).resolve() if not path.is_absolute() else path.resolve()
        require(path.is_relative_to(self.root), 'dependency outside project root')
        return path

    @staticmethod
    def signature(path):
        s = path.stat()
        require(path.is_file(), 'dependency is not a file')
        return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)

    def key(self, path):
        return os.path.normcase(str(path))

    def ensure_open(self):
        require(not self.finished, 'preflight transaction already sealed')

    def fresh_hash(self, path):
        before = self.signature(path)
        h = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024*1024), b''):
                h.update(block)
                self.hash_bytes += len(block)
        after = self.signature(path)
        require(before == after, 'dependency changed while hashing: '+str(path))
        self.hash_reads += 1
        return h.hexdigest(), after

    def sha(self, value):
        self.ensure_open(); self.requests += 1
        path = self.path(value); key = self.key(path)
        if key not in self.records:
            digest, signature = self.fresh_hash(path)
            self.records[key] = dict(path=path, digest=digest, signature=signature)
        row = self.records[key]
        require(self.signature(path) == row['signature'], 'dependency mutation within preflight: '+str(path))
        return row['digest']

    def verify(self, mapping):
        self.ensure_open()
        for name, digest in mapping.items():
            path = self.path(name); key = self.key(path)
            require(isinstance(digest, str) and len(digest) == 64,
                    'valid SHA256 expectation required: '+str(name))
            if key in self.expected:
                require(self.expected[key] == digest, 'conflicting expected hashes: '+str(name))
            self.expected[key] = digest
            require(self.sha(path) == digest, 'source/evidence drift: '+str(name))

    def read(self, value):
        path = self.path(value); key = self.key(path); digest = self.sha(path)
        if key not in self.json_cache:
            payload = path.read_bytes()
            require(hashlib.sha256(payload).hexdigest() == digest,
                    'JSON changed between hash and read: '+str(path))
            require(self.signature(path) == self.records[key]['signature'], 'JSON mutation during read')
            self.json_cache[key] = json.loads(payload.decode('utf-8'))
        return copy.deepcopy(self.json_cache[key])

    def load(self, value):
        path = self.path(value); key = self.key(path); digest = self.sha(path)
        if key not in self.load_cache:
            pair = original.load(path)
            require(self.signature(path) == self.records[key]['signature'], 'checkpoint changed during load')
            self.load_cache[key] = pair
        pa, state = self.load_cache[key]
        # Match the original repeated load semantics: no mutable ParticleArray
        # or metadata is shared with a different administrative caller.
        return pa.extract_particles(np.arange(len(pa.x), dtype=np.uint32)), copy.deepcopy(state)

    def checked_receipt(self, folder, statuses, frozen=None):
        path = self.path(folder); key = self.key(path)
        if key not in self.cpu_receipts:
            self.cpu_receipts[key] = self._checked_receipt_uncached(path, statuses, frozen)
        meta = self.cpu_receipts[key]
        require(meta['status'] in statuses, 'cached CPU receipt status not permitted by caller')
        if frozen is not None:
            require(meta['source_sha256'] == frozen, 'cached CPU source binding changed')
        self.sha(path/'manifest.json'); self.sha(path/'supervision.json')
        return copy.deepcopy(meta)

    def gpu_receipt(self, folder, statuses):
        path = self.path(folder); key = self.key(path)
        if key not in self.gpu_receipts:
            self.gpu_receipts[key] = self._gpu_receipt_uncached(path, statuses)
        meta = self.gpu_receipts[key]
        require(meta['status'] in statuses, 'cached GPU receipt status not permitted by caller')
        self.sha(path/'manifest.json'); self.sha(path/'supervision.json')
        return copy.deepcopy(meta)

    def finish(self):
        self.ensure_open()
        # A complete second byte-level pass detects even same-size/content
        # rewrites with restored mtime. Metadata alone never authorizes a run.
        for row in self.records.values():
            digest, signature = self.fresh_hash(row['path'])
            require(digest == row['digest'] and signature == row['signature'],
                    'dependency changed before final seal: '+str(row['path']))
        self.finished = True
        return dict(status='PASS_FRESH_SHA256_TRANSACTION_AND_FINAL_REHASH',
            unique_files=len(self.records), expected_paths=len(self.expected),
            hash_requests=self.requests, byte_hash_reads=self.hash_reads,
            hashed_bytes=self.hash_bytes, parsed_json_files=len(self.json_cache),
            decoded_checkpoints=len(self.load_cache), validated_cpu_receipts=len(self.cpu_receipts),
            validated_gpu_receipts=len(self.gpu_receipts), wall_seconds=time.perf_counter()-self.started,
            persistent_cache=False, stat_only_trust=False, original_module_globals_modified=False)

    # STATIC_ADMINISTRATIVE_METHODS_FOLLOW


    def _checked_receipt_uncached(self, folder, statuses, frozen=None):
        meta, supervision = (self.read(folder / 'manifest.json'), self.read(folder / 'supervision.json'))
        require(meta['status'] in statuses and supervision['exit_code'] == 0 and supervision['source_inputs_still_match'] and (supervision['manifest_sha256'] == self.sha(folder / 'manifest.json')) and (supervision['log_sha256'] == self.sha(folder / 'run.log')), 'stage receipt absent, failed or changed')
        self.verify(meta['source_sha256'])
        self.verify(meta['artifacts_sha256'])
        if frozen is not None:
            require(meta['source_sha256'] == frozen, 'new implementation requires a new independent cost stage')
        require(meta['classification'] == model.CLASSIFICATION and (not meta['accepted']) and (not meta['restart_allowed']) and (meta['accepted_time_advanced_s'] == 0.0), 'only unchanged deposition diagnostic states may continue')
        return meta

    def bind_receipt(self, folder, meta, references):
        references.update(meta['artifacts_sha256'])
        for name in ('manifest.json', 'supervision.json', 'run.log'):
            path = folder / name
            references[rel(path)] = self.sha(path)

    def continuation_parent(self, cost_folder, parent_folder, frozen):
        cost = self.checked_receipt(cost_folder, {COST_PASS}, frozen)
        references = {}
        self.bind_receipt(cost_folder, cost, references)
        initial = ROOT / cost['initial_checkpoint']
        require(self.sha(initial) == cost['initial_sha256'], 'independent initial checkpoint changed')
        target = initial
        parent_manifest = None
        visited = set()
        current = parent_folder
        while current is not None:
            current = current.resolve()
            require(current not in visited, 'cycle in diagnostic continuation history')
            visited.add(current)
            meta = self.checked_receipt(current, {STAGE_PASS, SHOULDER_PASS, PREFIX_PASS}, frozen)
            require(meta['screened_continuation_authorized'] and meta['screened_history'] and (meta['initial_sha256'] == cost['initial_sha256']) and (meta['initial_checkpoint'] == cost['initial_checkpoint']), 'stage origin or screening authorization differs')
            require(meta['screened_continuation_checkpoint'] == meta['screened_history'][-1]['checkpoint'] and meta['screened_continuation_sha256'] == meta['screened_history'][-1]['sha256'], 'only the last screened fine state may continue')
            require(len(meta['screened_history']) == 2 * len(meta['pair_screens']) and all((row['passed'] for row in meta['pair_screens'])), 'all selected fine pairs must pass')
            for row in meta['screened_history']:
                require(self.sha(ROOT / row['checkpoint']) == row['sha256'] and self.sha(ROOT / row['observer_path']) == row['observer_sha256'] and (self.sha(ROOT / row['screen_path']) == row['screen_sha256']), 'screened history evidence changed')
            for offset in range(0, len(meta['screened_history']), 2):
                first, last = meta['screened_history'][offset:offset + 2]
                require(first['half_index'] == 1 and last['half_index'] == 2 and (first['macro_index'] == last['macro_index']) and (first['screen_path'] == last['screen_path']), 'screened pair row identity mismatch')
                actual_screen = self.read(ROOT / first['screen_path'])
                require(actual_screen == meta['pair_screens'][offset // 2] and actual_screen['passed'] and actual_screen['global_screen']['pass'] and actual_screen['local_screen']['passed'] and actual_screen['crossing_ids_match'] and (actual_screen['fine_checkpoint_sha256'] == [first['sha256'], last['sha256']]) and (actual_screen['time_s'] == last['time_s']), 'PASS screen does not bind these exact fine candidates')
            if parent_manifest is None:
                parent_manifest = meta
                target = ROOT / meta['screened_continuation_checkpoint']
            self.bind_receipt(current, meta, references)
            prior = meta.get('parent_stage_manifest')
            if prior is None:
                require(meta['parent_checkpoint_sha256'] == cost['initial_sha256'], 'first stage must start at independent initial state')
                current = None
            else:
                prior_path = ROOT / prior
                require(self.sha(prior_path) == meta['parent_stage_manifest_sha256'], 'parent stage manifest changed')
                prior_meta = self.read(prior_path)
                require(meta['parent_checkpoint_sha256'] == prior_meta['screened_continuation_sha256'], 'stage chain skipped or mixed parents')
                current = prior_path.parent
        pa, state = self.load(target)
        model.validate_state(state)
        require(state['classification'] == model.CLASSIFICATION and state['diagnostic_only'] and (not state['accepted']) and (not state['restart_allowed']), 'forbidden formal/other-model parent')
        if parent_manifest is not None:
            require(state['time_s'] == parent_manifest['screened_end_time_s'] and state['step'] == parent_manifest['screened_end_step'], 'parent checkpoint endpoint differs')
        return (cost, target, pa, state, references, parent_manifest)

    def _gpu_receipt_uncached(self, folder, statuses):
        meta, supervision = (self.read(folder / 'manifest.json'), self.read(folder / 'supervision.json'))
        require(meta['status'] in statuses and supervision['exit_code'] == 0 and supervision['source_inputs_still_match'] and (supervision['manifest_sha256'] == self.sha(folder / 'manifest.json')) and (supervision['log_sha256'] == self.sha(folder / 'run.log')), 'invalid GPU receipt')
        self.verify(meta['source_sha256'])
        self.verify(meta['reference_sha256'])
        self.verify(meta['artifacts_sha256'])
        require(not meta['accepted'] and (not meta['restart_allowed']) and meta['diagnostic_only'] and (meta['classification'] == model.CLASSIFICATION) and (meta['accepted_time_advanced_s'] == 0.0), 'diagnostic GPU states only')
        if meta['status'] in {REPLAY_PASS, CONTINUE_PASS}:
            self.validate_gpu_screened_receipt(meta, folder)
        return meta

    def validate_gpu_screened_receipt(self, meta, folder):
        history, screen = (meta['screened_history'], meta['pair_screen'])
        require(meta['screened_continuation_authorized'] and len(history) == 2 and screen['passed'] and screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'], 'complete passing GPU fine pair required')
        require(self.read(folder / 'pair_screen.json') == screen, 'GPU manifest and pair screen differ')
        require(screen['fine_checkpoint_sha256'] == [r['sha256'] for r in history] and meta['screened_continuation_checkpoint'] == history[-1]['checkpoint'] and (meta['screened_continuation_sha256'] == history[-1]['sha256']), 'GPU continuation must bind exact last fine state')
        parent_pa, parent_state = self.load(ROOT / meta['parent_checkpoint'])
        require(self.sha(ROOT / meta['parent_checkpoint']) == meta['parent_sha256'] and cpu.parent_fingerprint(parent_pa, parent_state) == meta['parent_fingerprint'], 'GPU full parent identity changed')
        prefix_n = len(parent_state['machine_history'])
        for half, (role, row) in enumerate(zip(('fine_1', 'fine_2'), history), 1):
            require((ROOT / row['checkpoint']).resolve() == (folder / role / 'NOT_RESTARTABLE_candidate.npz').resolve() and self.sha(ROOT / row['checkpoint']) == row['sha256'], 'GPU fine checkpoint path/hash mismatch')
            _, state = self.load(ROOT / row['checkpoint'])
            model.validate_state(state)
            require(state['diagnostic_only'] and (not state['accepted']) and (not state['restart_allowed']) and (state['time_s'] == row['time_s']), 'GPU fine metadata mismatch')
            require(abs(state['time_s'] - meta['parent_time_s'] - half * meta['dt_s'] / 2) <= 1e-13 and state['step'] == parent_state['step'] + half, 'GPU half-step time/history discontinuity')
            require(state['machine_history'][:prefix_n] == parent_state['machine_history'] and state['interval_energy']['records'][:prefix_n] == parent_state['interval_energy']['records'], 'GPU parent history prefix differs')
        require(state['time_s'] == meta['screened_end_time_s'] and state['step'] == meta['screened_end_step'] and (abs(state['time_s'] - meta['parent_time_s'] - meta['dt_s']) <= 1e-13), 'GPU interval endpoint mismatch')
        require(self.sha(folder / 'coarse/NOT_RESTARTABLE_candidate.npz') == screen['coarse_checkpoint_sha256'], 'GPU coarse screen binding mismatch')

    def freeze(self, cpu_meta):
        mapping = dict(cpu_meta['source_sha256'])
        components_path = BASE / 'components_v1/manifest.json'
        components = self.read(components_path)
        require(components['status'] == 'PASS_GPU_COMPONENTS' and all((r['exit_code'] == 0 for r in components['checks'])), 'executed GPU components required')
        self.verify(components['source_sha256'])
        self.verify(components['artifacts_sha256'])
        mapping.update(components['source_sha256'])
        mapping.update(components['artifacts_sha256'])
        files = [CONTRACT, ORIGINAL_RUNNER, components_path]
        files += list((ROOT / 'scripts').glob('*afsd_gpu*20261007.py'))
        files += list(gpu.OVERLAY.rglob('*.pyd')) + list(gpu.OVERLAY.rglob('METADATA'))
        files += list(gpu.OVERLAY.rglob('*.py'))
        for path in files:
            mapping[rel(path)] = self.sha(path)
        self.verify(mapping)
        require(self.read(BASE / 'components_v1/linalg.json')['status'] == 'PASS_CUDA_FP64_COMPONENTS', 'GPU linalg checks required')
        require(self.read(BASE / 'components_v1/thermal.json')['status'] == 'PASS' and self.read(BASE / 'components_v1/thermal.json')['gpu_requested'], 'executed GPU thermal checks required')
        adaptation_audit()
        return mapping

    def references(self, folder, meta, mapping):
        mapping.update(meta['artifacts_sha256'])
        for name in ('manifest.json', 'supervision.json', 'run.log'):
            mapping[rel(folder / name)] = self.sha(folder / name)

    def preflight(self, args):
        baseline = (ROOT / args.baseline).resolve()
        require(baseline == (CPU_BASE / 'stage_013').resolve(), 'this registered replay is specifically CPU stage013')
        original = self.checked_receipt(baseline, {cpu.STAGE_PASS, cpu.PREFIX_PASS, cpu.SHOULDER_PASS})
        require(len(original['pair_screens']) == 1 and original['pair_screens'][0]['refinement_level'] == 0 and (original['screened_continuation_sha256'] == 'c3a9e4e5b194de04a5528b06d7fce6815c97a4f8b9d823aa53206dc773ff73e2'), 'registered exact CPU one-pair reference changed')
        cost_folder = ROOT / original['cost_directory']
        _, cpu_target, latest_pa, latest_state, refs, _ = self.continuation_parent(cost_folder, baseline, original['source_sha256'])
        self.references(baseline, original, refs)
        frozen = self.freeze(original)
        dt = original['pair_screens'][0]['dt_s']
        replay_parent = ROOT / original['parent_checkpoint']
        protocol = dict(schema='EXPLICIT_GPU_FP64_BRIDGE_20261007', mode=args.mode, source_sha256=frozen, reference_sha256=refs, classification=model.CLASSIFICATION, cpu_baseline=rel(baseline), cpu_stage_wall_s=original['wall_seconds'], cpu_candidate_wall_s=original['maximum_candidate_wall_s'], dt_s=dt, parent_checkpoint=rel(replay_parent), parent_sha256=self.sha(replay_parent), cpu_latest_checkpoint=rel(cpu_target), cpu_latest_sha256=self.sha(cpu_target), particle_count=len(latest_pa.x), free_dof=int(latest_state['free'].sum()), accepted=False, restart_allowed=False, diagnostic_only=True, accepted_time_advanced_s=0.0, internal_limit_s=540.0, external_limit_s=600.0, forecast_limit_s=450.0, compute_protocol='CPU_PYSPH_WITH_EXPLICIT_CUDA_FP64_DENSE_SOLVES_V1')
        if args.mode == 'benchmark':
            protocol['forecast_s'] = 90.0
        else:
            require(args.cost, '--cost required')
            cost_dir = (ROOT / args.cost).resolve()
            cost = self.gpu_receipt(cost_dir, {BENCH_PASS})
            require(cost['source_sha256'] == frozen and cost['cpu_baseline'] == rel(baseline), 'GPU benchmark binding changed')
            self.references(cost_dir, cost, refs)
            protocol['cost_directory'] = rel(cost_dir)
            if args.mode == 'replay':
                protocol['forecast_s'] = 30.0 + 1.5 * original['wall_seconds']
                protocol['candidate_forecast_s'] = original['maximum_candidate_wall_s']
            else:
                require(args.validation, '--validation required')
                valid_dir = (ROOT / args.validation).resolve()
                valid = self.gpu_receipt(valid_dir, {REPLAY_PASS})
                require(valid['source_sha256'] == frozen and valid['cpu_baseline'] == rel(baseline), 'GPU replay binding changed')
                self.references(valid_dir, valid, refs)
                chosen = valid
                if args.parent:
                    parent_dir = (ROOT / args.parent).resolve()
                    visited = set()
                    while True:
                        require(parent_dir not in visited, 'GPU continuation cycle')
                        visited.add(parent_dir)
                        prior = self.gpu_receipt(parent_dir, {CONTINUE_PASS})
                        require(prior['source_sha256'] == frozen and prior['validation_directory'] == rel(valid_dir), 'mixed GPU continuation')
                        self.references(parent_dir, prior, refs)
                        if chosen is valid:
                            chosen = prior
                        previous = prior.get('parent_gpu_directory')
                        earlier = self.gpu_receipt(ROOT / previous, {CONTINUE_PASS}) if previous else valid
                        require(prior['parent_sha256'] == earlier['screened_continuation_sha256'], 'GPU parent chain mismatch')
                        if not previous:
                            break
                        parent_dir = ROOT / previous
                protocol.update(parent_checkpoint=chosen['screened_continuation_checkpoint'], parent_sha256=chosen['screened_continuation_sha256'], validation_directory=rel(valid_dir), parent_gpu_directory=rel(ROOT / args.parent) if args.parent else None, dt_s=args.dt if args.dt is not None else chosen['dt_s'], candidate_forecast_s=chosen['maximum_candidate_wall_s'])
                require(0 < protocol['dt_s'] <= 0.005, 'bounded positive GPU macro interval required')
                protocol['forecast_s'] = 30.0 + 4.5 * protocol['candidate_forecast_s']
        require(protocol['forecast_s'] <= 450.0, 'GPU stage measured forecast exceeds450s')
        self.verify(refs)
        return protocol


STATIC_METHODS = (
    ('checked_receipt', '_checked_receipt_uncached', cpu),
    ('bind_receipt', 'bind_receipt', cpu),
    ('continuation_parent', 'continuation_parent', cpu),
    ('gpu_receipt', '_gpu_receipt_uncached', original),
    ('validate_gpu_screened_receipt', 'validate_gpu_screened_receipt', original),
    ('freeze', 'freeze', original),
    ('references', 'references', original),
    ('preflight', 'preflight', original),
)


def static_adaptation_audit(new_tree=None):
    """Every original administrative assertion survives the named IO changes."""
    tree = ast.parse(Path(__file__).read_text(encoding='utf-8')) if new_tree is None else copy.deepcopy(new_tree)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'CachedPreflight')
    methods = {n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)}
    relative_file = ast.dump(ast.parse('Path(__file__).resolve()', mode='eval').body, include_attributes=False)
    class NormalizeReference(ast.NodeTransformer):
        def visit_Call(self, node):
            if ast.dump(node, include_attributes=False) == relative_file:
                return ast.Name(id='ORIGINAL_RUNNER', ctx=ast.Load())
            f = node.func
            names = {'read', 'load', 'verify', 'checked_receipt', 'bind_receipt',
                'continuation_parent', 'gpu_receipt', 'validate_gpu_screened_receipt', 'freeze', 'references'}
            name = f.id if isinstance(f, ast.Name) and f.id in names else None
            if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name):
                if f.value.id == 'core' and f.attr == 'sha':
                    name = 'sha'
                if f.value.id == 'cpu' and f.attr in {'checked_receipt', 'continuation_parent'}:
                    name = f.attr
            if name:
                node.func = ast.Attribute(value=ast.Name(id='self', ctx=ast.Load()), attr=name, ctx=ast.Load())
            return self.generic_visit(node)
    for old_name, new_name, provider in STATIC_METHODS:
        reference = ast.parse(textwrap.dedent(inspect.getsource(getattr(provider, old_name)))).body[0]
        reference = NormalizeReference().visit(reference)
        reference.name = new_name; reference.args.args.insert(0, ast.arg(arg='self'))
        require(new_name in methods and ast.dump(reference, include_attributes=False) ==
            ast.dump(methods[new_name], include_attributes=False), 'administrative gate changed: '+new_name)
    require(model.validate_state is cpu.model.validate_state, 'model validation binding changed')
    return dict(status='PASS_STATIC_ADMINISTRATIVE_IO_ADAPTATION_ONLY', methods=len(STATIC_METHODS),
        all_original_assertions_retained=True, physics_worker=str(ORIGINAL_RUNNER),
        no_physical_source_modified=True, no_original_module_global_rebinding=True)


def bind_administration(transaction, protocol, checks_path):
    checks_path = transaction.path(checks_path)
    checks = transaction.read(checks_path)
    require(checks['status'] == 'PASS_GPU_PREFLIGHT_CACHE_COMPONENTS' and checks['checks']
        and all(row['exit_code'] == 0 for row in checks['checks']), 'executed administrative component checks required')
    transaction.verify(checks['source_sha256']); transaction.verify(checks['artifacts_sha256'])
    required = {rel(Path(__file__)), rel(ADMIN_CONTRACT), 'tests/review/test_gpu_preflight_cache_20261008.py'}
    require(required.issubset(checks['source_sha256']), 'administrative test receipt misses required sources')
    admin = dict(checks['source_sha256'])
    for mapping in (checks['artifacts_sha256'], {rel(checks_path): transaction.sha(checks_path)}):
        for key, value in mapping.items():
            require(key not in admin or admin[key] == value, 'administrative source/artifact hash conflict')
            admin[key] = value
    for key, value in admin.items():
        require(key not in protocol['reference_sha256'] or protocol['reference_sha256'][key] == value,
                'administrative and physical reference hash conflict')
        require(key not in protocol['source_sha256'] or protocol['source_sha256'][key] == value,
                'administrative and physical source hash conflict')
        protocol['reference_sha256'][key] = value
    protocol.update(administrative_source_sha256=admin,
        administrative_protocol='PER_LAUNCH_FRESH_BYTE_HASH_DEDUPLICATION_WITH_FINAL_REHASH_V1',
        administrative_checks=rel(checks_path), administrative_adaptation_audit=static_adaptation_audit())


def create_protocol(args, bind_admin=True):
    transaction = CachedPreflight()
    static_adaptation_audit()
    protocol = transaction.preflight(args)
    if bind_admin:
        bind_administration(transaction, protocol, args.admin_checks)
    stats = transaction.finish()
    return protocol, stats


def launch_original_worker(args, protocol, statistics):
    output = (ROOT/args.output).resolve()
    require(output.is_relative_to(BASE) and output != BASE and not output.exists(),
            'new separate GPU output required')
    output.mkdir(parents=True)
    protocol['administrative_preflight_statistics'] = statistics
    original.shared.write_new(output/'protocol.json', protocol)
    command = [sys.executable, '-B', str(ORIGINAL_RUNNER), '--mode', args.mode,
               '--output', rel(output), '--worker']
    env = dict(os.environ, **{key: '1' for key in original.shared.THREADS}); env['PYTHONPATH'] = ''
    began = time.perf_counter(); code = 124
    with (output/'run.log').open('x', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=ROOT, env=env, stdout=stream,
                stderr=subprocess.STDOUT, timeout=600.).returncode
        except subprocess.TimeoutExpired:
            stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    error = None; post_statistics = None
    try:
        # A new transaction after the child performs fresh byte verification.
        # Nothing from the earlier preflight can authorize an altered source.
        after = CachedPreflight()
        after.verify(protocol['source_sha256']); after.verify(protocol['reference_sha256'])
        post_statistics = after.finish()
    except Exception as exc:
        code = 2; error = repr(exc)
    manifest = output/'manifest.json'
    original.shared.write_new(output/'supervision.json', dict(exit_code=code, command=command,
        wall_seconds=time.perf_counter()-began, source_inputs_still_match=error is None, source_error=error,
        manifest_sha256=core.sha(manifest) if manifest.exists() else None,
        log_sha256=core.sha(output/'run.log'), accepted=False, restart_allowed=False,
        administrative_launcher=rel(Path(__file__)), administrative_source_sha256=protocol['administrative_source_sha256'],
        administrative_preflight_statistics=statistics, administrative_postflight_statistics=post_statistics))
    print(json.dumps(dict(exit_code=code, output=rel(output),
        status=original.read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    return code


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--mode', default='continue', choices=('continue',))
    result.add_argument('--output')
    result.add_argument('--baseline', default=rel(CPU_BASE/'stage_013'))
    result.add_argument('--cost', default=rel(BASE/'benchmark_v1'))
    result.add_argument('--validation', default=rel(BASE/'replay_v1'))
    result.add_argument('--parent'); result.add_argument('--dt', type=float)
    result.add_argument('--admin-checks', type=Path, default=DEFAULT_CHECKS)
    result.add_argument('--preflight-only', action='store_true')
    result.add_argument('--compare-protocol', type=Path)
    result.add_argument('--comparison-output', type=Path)
    return result


def main():
    args = parser().parse_args()
    if args.compare_protocol:
        require(args.preflight_only, 'protocol comparison is read-only preflight only')
        path = (ROOT/args.compare_protocol).resolve(); saved = original.read(path)
        args.baseline, args.cost, args.validation = saved['cpu_baseline'], saved['cost_directory'], saved['validation_directory']
        args.parent, args.dt = saved.get('parent_gpu_directory'), saved['dt_s']
        protocol, statistics = create_protocol(args, bind_admin=False)
        require(protocol == saved, 'cached preflight differs from saved original protocol')
        report = dict(status='PASS_IDENTICAL_SAVED_PHYSICAL_PROTOCOL_NO_SOLVER',
            reference_protocol=rel(path), reference_sha256=core.sha(path), statistics=statistics,
            physical_time_advanced_s=0., physics_executed=False, gpu_executed=False)
        if args.comparison_output:
            destination = (ROOT/args.comparison_output).resolve()
            require(destination.is_relative_to(ROOT) and not destination.exists(), 'new project-local comparison receipt required')
            destination.parent.mkdir(parents=True, exist_ok=True)
            original.shared.write_new(destination, report)
        print(json.dumps(report)); return 0
    require(args.output, '--output required')
    output = (ROOT/args.output).resolve()
    require(output.is_relative_to(BASE) and output != BASE and not output.exists(), 'new GPU output required')
    protocol, statistics = create_protocol(args)
    if args.preflight_only:
        print(json.dumps(dict(status='PASS_CACHED_ADMINISTRATIVE_PREFLIGHT_ONLY',
            parent_checkpoint=protocol['parent_checkpoint'], dt_s=protocol['dt_s'], statistics=statistics)))
        return 0
    return launch_original_worker(args, protocol, statistics)


if __name__ == '__main__':
    raise SystemExit(main())
