"""Bounded DEBUG/APPROXIMATE host entry; production remains blocked."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import platform
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import threading
import numpy as np
from afsd_scaffold.output_contract import validate_manifest
from .config_adapter import load_settings
from .particles import create_geometry
from .output import vtk, checkpoint, restore, digest
from .solver import SolidSolver
from .geometry import tool_pose

PIN = '1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5'


def json_write(path, value):
    target = Path(path)
    temporary = target.with_suffix(target.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    os.replace(temporary, target)


def memory_bytes():
    if sys.platform == 'win32':
        import ctypes
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD)]+[(k, ctypes.c_size_t) for k in ('PeakWorkingSetSize', 'WorkingSetSize', 'QuotaPeakPagedPoolUsage', 'QuotaPagedPoolUsage', 'QuotaPeakNonPagedPoolUsage', 'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage')]
        value = Counters(); value.cb = ctypes.sizeof(value)
        kernel = ctypes.WinDLL('kernel32'); psapi = ctypes.WinDLL('psapi')
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(value), value.cb):
            raise OSError('GetProcessMemoryInfo failed')
        return int(value.PeakWorkingSetSize)
    import resource
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', default='configs/debug/starved.json')
    parser.add_argument('--output', required=True)
    parser.add_argument('--steps', type=int, default=200, help='additional steps, not total on restart')
    parser.add_argument('--output-every', type=int, default=20)
    parser.add_argument('--end-time', type=float, help='absolute SI terminal time; last dt shortened exactly')
    parser.add_argument('--wall-seconds', type=float, default=None)
    parser.add_argument('--restart', help='trusted project NPZ only (upstream uses pickle)')
    parser.add_argument('--contact-variant', choices=('legacy', 'event_aware_v1'), default='legacy',
                        help='explicit APPROXIMATE local event-aware kick; legacy is retained')
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[3]
    config, p, consumed = load_settings(args.case)
    classification = config['intended_run_classification']
    if args.end_time is not None and (not np.isfinite(args.end_time) or args.end_time <= 0):
        raise ValueError('positive finite end-time required')
    if args.steps <= 0 or args.output_every <= 0: raise ValueError('positive step counts required')
    wall_limit = p['wall_seconds'] if args.wall_seconds is None else args.wall_seconds
    if not (0 < wall_limit <= 600): raise ValueError('this runner caps trials at 600s; longer work needs explicit reviewed runner change')
    repository = root/'third_party/PySPH'
    commit = subprocess.check_output(['git', '-C', str(repository), 'rev-parse', 'HEAD'], text=True).strip()
    dirty = subprocess.check_output(['git', '-C', str(repository), 'status', '--porcelain', '--untracked-files=no'], text=True).strip()
    if commit != PIN or dirty: raise RuntimeError('upstream pin/cleanliness gate failed')
    source_paths = list((root/'src/afsd_backend').rglob('*.py'))+[root/'src/afsd_scaffold/config.py', root/'src/afsd_scaffold/vtk_legacy.py', root/'src/afsd_scaffold/output_contract.py', root/'scripts/run_afsd_debug.py']
    source_hashes = {str(f.relative_to(root)): digest(f) for f in sorted(source_paths)}
    fingerprint = hashlib.sha256(json.dumps(config['parameters'], sort_keys=True).encode()).hexdigest()
    if args.restart:
        pa, state = restore(args.restart)
        if state.get('configuration_fingerprint') != fingerprint or state.get('source_hashes') != source_hashes:
            raise RuntimeError('restart configuration/source mismatch; explicit migration required')
    else:
        pa = create_geometry(p); state = None
    if pa.get_number_of_particles() > 10000: raise ValueError('host DEBUG runner limited to 10000 particles')
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    config_path = output/'resolved_config.json'; json_write(config_path, config)
    json_write(output/'consumed_parameters.json', consumed)
    json_write(output/'source_hashes.json', source_hashes)
    solver = SolidSolver(pa, p, state, contact_variant=args.contact_variant)
    if args.end_time is not None and args.end_time <= solver.time:
        raise ValueError('end-time must exceed checkpoint time')
    log_file = (output/'run.log').open('w', encoding='utf-8')
    diagnostics_file = (output/'diagnostics.jsonl').open('w', encoding='utf-8')
    def log(message):
        print(message, flush=True); log_file.write(message+'\n'); log_file.flush()
    frames = []; fields = []
    def save_frame():
        nonlocal fields
        name = f'{classification}_state_{solver.count:06d}.vtk'
        fields = vtk(output/name, pa, p, solver.bulk, solver.time, classification=classification)
        position, angle = tool_pose(solver.time, p['traverse'], p['omega'], p['layer'])
        frames.append(dict(step=solver.count, time_s=solver.time, file=name, sha256=digest(output/name), fields=fields, rate_evaluation_time_s=solver.time-solver.last['dt_s']/2 if solver.count else None, unavailable_fields_reason='initial rate fields have not been evaluated' if not solver.count else None, tool_translation_m=position.tolist(), tool_rotation_rad=angle, tool_axis=[0, 0, 1]))
        checkpoint(output/'checkpoint.npz', pa, {**solver.state(), 'configuration_fingerprint': fingerprint, 'source_hashes': source_hashes})
        json_write(output/'series.json', {'classification': classification, 'frames': frames})
    log('command: '+subprocess.list2cmdline(sys.argv))
    log(f'backend=PySPH {commit}; host NumPy equations; no GPU, no JIT, no scaling; N={pa.get_number_of_particles()}')
    log(f'dt_initial_s={solver.stable_dt():.9g}; layer/dx={p["layer"]/p["dx"]}; wall_cap_s={wall_limit}')
    log('contact_variant='+args.contact_variant+'; production_authorized=false')
    save_frame()
    def make_manifest(status):
        manifest = dict(schema_version='0.1.0', run_id=output.name, run_classification=classification, physics_claims_allowed=True, created_utc=datetime.now(timezone.utc).isoformat(), status=status, backend=dict(name='PySPH + REVIEW_APPROX host adapter', repository='https://github.com/pypr/pysph', commit=commit, dirty=False, compiler=platform.python_compiler(), build_mode='pinned compiled ParticleArray/LinkedListNNPS + unmodified host LeapFrogIntegrator/Step; NumPy equations; CPU'), configuration=dict(case_id=config['case_id'], path='resolved_config.json', sha256=digest(config_path)), units=dict(calculation_system='SI', length='m', time='s', mass='kg', temperature='K', stress='Pa'), coordinate_frames=dict(stored_frame='LAB', definition='+x traverse, +z substrate normal, omega +z CCW; translating frame x-Xtool (no rotation); AS y<0, RS y>0', tool_pose_file='series.json'), particles=dict(file=last['file'], format='VTK_LEGACY_ASCII', sha256=last['sha256'], count=pa.get_number_of_particles(), persistent_id_field='particle_id', initial_position_field='initial_position_m', material_id_field='material_id', tracer_type_field='tracer_type', fields=fields), provenance=dict(command=subprocess.list2cmdline(sys.argv), logs=['run.log', 'diagnostics.jsonl', 'performance.json'], source_files=list(source_hashes), approximations=['docs/PARAMETER_SOURCES.md D01-D10', 'preheated finite bar, scaled lateral geometry, pre-bonded SPH interface, no startup/lift simulation', 'dimensionless b1; explicit rate floor; constant thermal properties; linear bulk EOS', 'finite analytic rigid faces, midpoint-work heat, memoryless sliding; incomplete rigid exterior contact', 'first-order constitutive/thermal split, corrected objective/work-conjugate SPH, no stabilization, no mass/time scaling']), checkpoint=dict(file='checkpoint.npz', sha256=digest(output/'checkpoint.npz'), step=last['step']), physics_scope='short thermomechanical implementation probe only; not deposition/mixing/experimental validation')
        manifest['provenance']['checkpoint_parent'] = str(Path(args.restart).resolve()) if args.restart else None
        contract = root/'logs/review/contract_v1_pretest.md'
        manifest['validation_contract'] = {'version': '1', 'path': str(contract.relative_to(root)), 'sha256': digest(contract), 'production_authorized': False}
        manifest['model_fidelity'] = 'APPROXIMATE; scientific review blocks production'
        manifest['provenance']['contact_variant'] = args.contact_variant
        if args.contact_variant == 'event_aware_v1':
            manifest['provenance']['approximations'][3] = 'finite analytic rigid faces; frozen-geometry/internal-force event-aware coupled kick, midpoint subinterval heat; incomplete rigid exterior contact'
        return manifest
    last = frames[-1]
    initial_manifest = make_manifest('INCOMPLETE')
    validate_manifest(initial_manifest)
    json_write(output/'manifest.json', initial_manifest)
    start = time.perf_counter(); step_times = []; dts = []
    status = 'COMPLETE'; error = None
    starting_count = solver.count; starting_time = solver.time
    def hard_timeout():
        json_write(output/'HARD_TIMEOUT.json', dict(status='INCOMPLETE', limit_seconds=wall_limit,
            classification='REJECTED_IN_PROGRESS', last_checkpoint_step=frames[-1]['step']))
        os._exit(124)
    watchdog = threading.Timer(wall_limit, hard_timeout)
    watchdog.daemon = True
    watchdog.start()
    try:
        for _ in range(args.steps):
            if args.end_time is not None and solver.time >= args.end_time:
                break
            if time.perf_counter()-start >= wall_limit:
                status = 'INCOMPLETE'; error = 'wall-clock cap reached at step boundary'; break
            remaining = args.end_time-solver.time if args.end_time is not None else None
            tick = time.perf_counter(); dt = solver.advance(max_dt=remaining)
            dts.append(dt); step_times.append(time.perf_counter()-tick)
            diag = solver.diagnostics(fail=False)
            diagnostics_file.write(json.dumps(diag, allow_nan=False)+'\n'); diagnostics_file.flush()
            if diag['alarms']: raise RuntimeError('; '.join(diag['alarms']))
            if solver.count % args.output_every == 0:
                save_frame()
                log(f'step={solver.count} t={solver.time:.9g}s energy_residual={diag["energy_residual_relative"]:.5g} ep_max={diag["ep_max"]:.5g}')
        if args.end_time is not None and solver.time < args.end_time and status == 'COMPLETE':
            status = 'INCOMPLETE'; error = 'step cap reached before requested physical end-time'
        if not frames or frames[-1]['step'] != solver.count: save_frame()
    except Exception as exception:
        status = 'FAILED'; error = traceback.format_exc(); log(error)
        # Last completed output checkpoint stays usable; no partial-step overwrite.
        checkpoint(output/'REJECTED_candidate.npz', pa, {**solver.state(),
            'classification':'REJECTED', 'restart_allowed':False, 'reason':error,
            'source_hashes':source_hashes, 'configuration_fingerprint':fingerprint,
            'stage':'partial_or_alarm_failed_candidate; last accepted endpoint is checkpoint.npz'})
        from .kick_diagnostics import clean
        json_write(output/'REJECTED_contact_debug.json', clean(dict(
            classification='REJECTED', restart_allowed=False,
            last_completed_contact_callback_trace=solver.last_contact_trace,
            # A domain_kick exception leaves the previous successful trace.
            # Treat this auxiliary trace conservatively; the separate local
            # failure context, its time and particle ID identify the candidate.
            trace_not_guaranteed_to_be_current_candidate=True,
            failed_local_context=getattr(exception, 'contact_context', None),
            last_accepted_checkpoint_step=frames[-1]['step'], error=error)))
    watchdog.cancel()
    elapsed = time.perf_counter()-start
    count = solver.count-starting_count
    physical = solver.time-starting_time
    report = dict(status=status, error=error, completed_steps=count, requested_steps=args.steps, wall_seconds=elapsed, wall_limit_seconds=wall_limit, physical_seconds=physical, particle_count=pa.get_number_of_particles(), peak_process_working_set_bytes=memory_bytes(), dt_min_s=min(dts) if dts else None, dt_max_s=max(dts) if dts else None, physics_step_median_seconds=float(np.median(step_times)) if step_times else None, end_to_end_steps_per_second=count/elapsed, physical_seconds_per_wall_second=physical/elapsed, final_checkpoint_step=frames[-1]['step'], classification=classification, quantitative_validation=False, requested_end_time_s=args.end_time)
    if dts:
        paper_dt = min(dts)*(.001/p['dx'])
        # Explicitly optimistic linear N scaling, not an observed production cost.
        report['paper_cost_projection'] = dict(particles=80312, nominal_spacing_m=.001, assumed_dt_s=paper_dt, model='optimistic linear N scaling of measured median host step; no claim for changed geometry/contact', wall_hours_per_physical_second=float(np.median(step_times))*80312/pa.get_number_of_particles()/paper_dt/3600, layer_spacings=1.5, warning='paper nominal spacing does not meet this debug layer>=3 gate; projection is not runnable validation')
    json_write(output/'performance.json', report)
    last = frames[-1]
    manifest = make_manifest(status)
    validate_manifest(manifest); json_write(output/'manifest.json', manifest)
    log(json.dumps(report)); diagnostics_file.close(); log_file.close()
    return 0 if status == 'COMPLETE' else (2 if status == 'INCOMPLETE' else 1)
