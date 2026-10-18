"""Independent Q5 saved-data audit; never imports the tested contact solver."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud
from pysph.solver.utils import load


def wall_reference(t):
    """Declared rigid boundary only, not a prescribed material trajectory."""
    x = -4e-5*(t/.05) if t <= .05 else -4e-5
    if t <= .05:
        y = 0.
    elif t <= .1:
        y = 4e-4*(t-.05)
    elif t <= .15:
        y = 2e-5-8e-4*(t-.1)
    else:
        y = -2e-5+4e-4*(t-.15)
    return np.array([x, y, 0.])


def same(actual, expected, tolerance, label):
    a, b = np.asarray(actual), np.asarray(expected)
    if a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=tolerance):
        raise ValueError(label)


def inspect(directory):
    path = Path(directory)
    r = json.loads((path/'report.json').read_text(encoding='utf-8'))
    cfg = json.loads((path/'resolved_config.json').read_text(encoding='utf-8'))
    required_sources = {'scripts/probe_process_friction.py', 'configs/benchmarks/process_friction_q5.json',
        'scripts/probe_process_contact.py', 'scripts/probe_process_plastic.py', 'scripts/probe_process_time.py',
        'scripts/probe_engineering_material.py', 'scripts/probe_engineering_block.py', 'scripts/probe_boundary_sbp.py',
        'src/afsd_backend/pysph/material.py', 'src/afsd_backend/pysph/mechanics.py'}
    if set(r['source_sha256']) != required_sources:
        raise ValueError('source manifest inventory mismatch')
    if not {'history.npz', 'checkpoint.npz', 'resolved_config.json', 'final.vtk'} <= set(r['artifact_sha256']):
        raise ValueError('artifact manifest inventory mismatch')
    for name, digest in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError('source changed: '+name)
    for name, digest in r['artifact_sha256'].items():
        if Path(name).name != name or hashlib.sha256((path/name).read_bytes()).hexdigest() != digest:
            raise ValueError('artifact changed or invalid path: '+name)
    if r['production_authorized'] is not False or cfg['production_authorized'] is not False:
        raise ValueError('Q5 cannot authorize production')
    if r['status'] not in ('DIRECTIONAL_FRICTION_PATH_COMPLETE', 'BLOCKED'):
        raise ValueError('unknown solver completion state')
    if r['increments_per_phase'] not in (4, 8) or cfg['increments_per_phase'] != r['increments_per_phase']:
        raise ValueError('undeclared Q5 requested grid')
    physics = dict(mu=.1, normal_displacement=-4e-5, tangent_amplitude=2e-5, phase_duration=.05, temperature=293.15)
    if any(cfg['parameters'][key]['value'] != value for key, value in physics.items()):
        raise ValueError('declared Q5 physical parameters changed')
    if not {'thermal', 'inertia', 'isotropic_two_tangent_friction', 'production_restart'} <= set(r['not_solved']):
        raise ValueError('missing unsolved physics declaration')
    with np.load(path/'history.npz', allow_pickle=False) as saved:
        h = {key: saved[key].copy() for key in saved.files}
    if not all(np.isfinite(a).all() for a in h.values()):
        raise ValueError('nonfinite accepted history')
    if any('temperature' in key or 'velocity' in key for key in h):
        raise ValueError('unsolved field published in history')
    t, rows = h['physical_time_s'], r['records']
    X, V, ids = h['reference_position_m'], h['reference_volume_m3'], h['particle_id']
    u, stress = h['displacement_m'], h['stress_Pa']
    if (len(t) != len(rows) or len(t) != r['accepted_steps']+1 or t[0] != 0 or
            t[-1] != r['accepted_time_s'] or np.any(np.diff(t) <= 0) or t[-1] > .2):
        raise ValueError('accepted frame/time inventory mismatch')
    grid = np.array([[a,b,c] for a in (0.,.005,.01) for b in (0.,.005,.01) for c in (0.,.005,.01)])
    if (X.shape != (27,3) or set(map(tuple, X)) != set(map(tuple, grid)) or ids.shape != (27,) or
            len(np.unique(ids)) != 27 or not np.equal(ids,np.floor(ids)).all()):
        raise ValueError('particle identity/reference geometry mismatch')
    nodes = np.flatnonzero(X[:,0] == .01)
    fixed = np.broadcast_to((X[:,0] == 0)[:,None], X.shape)
    if (not np.array_equal(h['support_components'], fixed) or np.any(u[:,fixed] != 0) or
            not np.array_equal(h['position_m'], X[None]+u) or r['permanent_free_dofs'] != 54 or
            r['particles'] != 27 or r['contact_rows'] != 9):
        raise ValueError('geometry/support/free-DOF mismatch')
    if (not np.array_equal(h['contact_particle_id'], ids[nodes]) or
            not np.array_equal(h['contact_feature_id'], np.ones(9, dtype=int))):
        raise ValueError('contact identity mismatch')
    expected_volume = .005**3*np.prod(np.where((X==0)|(X==.01), .5, 1.), axis=1)
    same(V, expected_volume, 1e-22, 'reference quadrature mismatch')
    same(h['mass_kg'], 2700*V, 1e-15, 'physical mass mismatch')
    if abs(h['mass_kg'].sum()/.0027-1) > 1e-12 or np.any(h['equivalent_plastic_strain'] != 0):
        raise ValueError('mass or declared elastic-only state mismatch')
    force, support = h['contact_force_N'], h['support_force_N']
    normal, tau, modes = h['normal_force_N'], h['tangent_force_N'], h['interval_modes']
    Q, dq = h['friction_dissipation_J'], h['interval_friction_J']
    if (normal.shape != (len(t),9) or tau.shape != normal.shape or modes.shape != normal.shape or
            Q.shape != normal.shape or dq.shape != normal.shape or not np.isin(modes,[-1,0,1]).all()):
        raise ValueError('contact history shape/mode mismatch')
    assembled = np.zeros_like(force)
    assembled[:,nodes,0] = -normal
    assembled[:,nodes,1] = tau
    same(force, assembled, 0, 'contact force assembly mismatch')
    if np.any(support[:,~fixed] != 0) or np.any(normal < -9e-5):
        raise ValueError('invalid support force or normal detachment')
    wall = np.array([wall_reference(time) for time in t])
    same(u[:,nodes,0], np.broadcast_to(wall[:,0,None], (len(t),9)), 1e-11, 'normal contact gap mismatch')
    if np.any(np.abs(tau)-.1*normal > 9e-5):
        raise ValueError('directional Coulomb bound violated')
    sliding = modes != 0
    same(tau[sliding], (-.1*modes*normal)[sliding], 1e-10, 'sliding traction sign/magnitude mismatch')
    for name in ('interval_friction_J', 'interval_relative_slip_m', 'tangent_impulse_N_s',
                 'relative_kinematic_rate_m_s', 'kinematic_interval_rate_m_s', 'friction_dissipation_J'):
        if np.any(h[name][0] != 0):
            raise ValueError('nonzero initial interval ledger')
    Ww = Wc = Wr = 0.
    reconstructed_Q = np.zeros(9)
    max_defect = max_force = max_torque = 0.
    for i, row in enumerate(rows):
        if row['step'] != i or row['physical_time_s'] != t[i]:
            raise ValueError('ledger/history frame mismatch')
        if row['mode_counts'] != {str(m):int((modes[i]==m).sum()) for m in (-1,0,1)}:
            raise ValueError('mode ledger mismatch')
        same(row['wall_displacement_m'], wall[i], 1e-19, 'wall program ledger mismatch')
        same(row['wall_force_on_material_N'], force[i].sum(axis=0), 1e-9, 'wall reaction ledger mismatch')
        max_force = max(max_force, float(np.linalg.norm((force[i]+support[i]).sum(axis=0))))
        max_torque = max(max_torque, float(np.linalg.norm(np.cross(X,force[i]+support[i]).sum(axis=0))))
        if i:
            dt = t[i]-t[i-1]
            du, dw = u[i]-u[i-1], wall[i]-wall[i-1]
            mid = .5*(force[i]+force[i-1])
            slip = du[nodes,1]-dw[1]
            impulse = mid[nodes,1]*dt
            relative_rate = slip/dt
            raw = -mid[nodes,1]*slip
            if np.any(raw < -1e-10):
                raise ValueError('raw negative per-contact dissipation')
            if (np.any(np.abs(slip[modes[i]==0]) > 1e-11) or
                    np.any(modes[i,sliding[i]]*slip[sliding[i]] < -1e-11)):
                raise ValueError('stick/slip direction violation')
            same(h['interval_relative_slip_m'][i], slip, 1e-18, 'relative slip mismatch')
            same(h['tangent_impulse_N_s'][i], impulse, 1e-12, 'tangential impulse mismatch')
            # Equivalent wall formulas can differ in their last displacement bit;
            # dividing by a tiny event interval must not turn that into a false failure.
            same(h['relative_kinematic_rate_m_s'][i], relative_rate, 1e-15+1e-18/dt, 'relative rate mismatch')
            same(h['kinematic_interval_rate_m_s'][i], du/dt, 1e-15, 'kinematic rate mismatch')
            same(dq[i], raw, 1e-12, 'per-contact friction work mismatch')
            same(dq[i], -impulse*relative_rate, 1e-12, 'impulse-relative-work mismatch')
            wc = float(np.sum(mid*du))
            ww = float(mid.sum(axis=0)@dw)
            wr = float(np.sum(.5*(support[i]+support[i-1])*du))
            if abs(ww-wc-raw.sum()) > 1e-6:
                raise ValueError('interval wall/material/friction identity mismatch')
            same(row['interval_wall_work_J'], ww, 1e-12, 'wall interval ledger mismatch')
            same(row['interval_contact_material_work_J'], wc, 1e-12, 'material interval ledger mismatch')
            Ww += ww; Wc += wc; Wr += wr; reconstructed_Q += raw
        same(Q[i], reconstructed_Q, 1e-12, 'per-contact accumulated dissipation mismatch')
        same(row['friction_dissipation_J'], reconstructed_Q.sum(), 1e-12, 'summed friction ledger mismatch')
        same(row['interval_friction_J'], dq[i].sum(), 1e-12, 'interval friction sum mismatch')
        same(row['minimum_raw_contact_dissipation_J'], dq[i].min(), 1e-12, 'minimum raw friction mismatch')
        same(row['wall_work_J'], Ww, 1e-12, 'cumulative wall ledger mismatch')
        same(row['contact_material_work_J'], Wc, 1e-12, 'cumulative material ledger mismatch')
        same(row['support_work_J'], Wr, 1e-12, 'support work mismatch')
        defect = Ww-Wc-reconstructed_Q.sum()
        same(row['friction_work_identity_J'], defect, 1e-12, 'reported contact identity mismatch')
        max_defect = max(max_defect, abs(defect))
        if abs(defect) > 1e-6 or abs(Wc+Wr-row['stress_work_J']) > 1e-6:
            raise ValueError('cumulative work mismatch')
        trace = np.trace(stress[i], axis1=1, axis2=2)
        dev = stress[i]-trace[:,None,None]*np.eye(3)/3
        Es = float(np.sum(V*(np.sum(dev*dev,axis=(1,2))/(4*35e9)+trace**2/(18*(70e9/3)))))
        same(row['elastic_diagnostic_J'], Es, 1e-10, 'independent elastic diagnostic mismatch')
        if row['plastic_work_J'] != 0 or row['maximum_ep'] != 0:
            raise ValueError('elastic-only ledger has plasticity')
    event_steps = [event['step'] for event in r['events']]
    if event_steps != [i for i,row in enumerate(rows) if row['event_kind']=='FRICTION_LIMIT']:
        raise ValueError('event/ledger inventory mismatch')
    for event in r['events']:
        i, bracket = event['step'], event['bracket_s']
        if (event['kind'] != 'FRICTION_LIMIT' or not 0<i<len(t) or event['physical_time_s'] != t[i] or
                len(bracket) != 2 or bracket[0] != t[i] or not 0<=bracket[1]-bracket[0]<=1e-12):
            raise ValueError('accepted friction event mismatch')
        same(event['ending_interval_modes'], modes[i], 0, 'event mode identity mismatch')
        same(event['cone_margin_N'], .1*normal[i]-np.abs(tau[i]), 1e-10, 'event cone margin mismatch')
    vtk = read_point_cloud(path/'final.vtk')
    saved = load(str(path/'checkpoint.npz')); pa = saved['arrays']['solid']; state = saved['solver_data']
    if (state['restart_allowed'] is not False or state['production_authorized'] is not False or
            state['step'] != r['accepted_steps'] or state['physical_time_s'] != t[-1]):
        raise ValueError('diagnostic checkpoint metadata mismatch')
    for saved_ids in (vtk['scalars']['particle_id'], pa.particle_id):
        same(saved_ids, ids, 0, 'accepted particle ID mismatch')
    for tags in (vtk['scalars']['tracer_type'], pa.tracer_type):
        same(tags, h['tracer_type'], 0, 'accepted passive tag mismatch')
    same(vtk['points'], h['position_m'][-1], 0, 'VTK accepted position mismatch')
    same(np.column_stack((pa.x,pa.y,pa.z)), h['position_m'][-1], 0, 'checkpoint position mismatch')
    same(pa.ep, h['equivalent_plastic_strain'], 0, 'checkpoint plastic strain mismatch')
    same(pa.temperature, np.full(27,293.15), 0, 'prescribed temperature input mismatch')
    for a,b,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        same(vtk['scalars']['stress_'+name+'_Pa'], stress[-1,:,a,b], 0, 'VTK stress mismatch')
    if any('temperature' in key or 'velocity' in key for key in vtk['scalars']|vtk['vectors']):
        raise ValueError('unsolved field published in VTK')
    complete = r['status']=='DIRECTIONAL_FRICTION_PATH_COMPLETE'
    if complete:
        requested = np.linspace(0.,.2,4*r['increments_per_phase']+1)
        if t[-1] != .2 or any(not np.any(np.abs(t-time)<=1e-14) for time in requested):
            raise ValueError('complete path lacks terminal/requested state')
    else:
        rejected = json.loads((path/'REJECTED_manifest.json').read_text(encoding='utf-8'))
        if (rejected['classification'] != 'REJECTED' or rejected['restart_allowed'] is not False or
                rejected['last_accepted_step'] != r['accepted_steps'] or not r['error']):
            raise ValueError('BLOCKED/rejected state metadata mismatch')
    return dict(status='PASS_COMPLETE_READBACK' if complete else 'PASS_PARTIAL_READBACK',
        completion_status=r['status'], production_authorized=False, accepted_steps=r['accepted_steps'],
        accepted_time_s=float(t[-1]), events=r['events'], independent_wall_work_J=Ww,
        independent_material_contact_work_J=Wc, independent_friction_dissipation_J=reconstructed_Q.tolist(),
        maximum_cumulative_contact_work_defect_J=max_defect, minimum_raw_contact_increment_J=float(dq.min()),
        maximum_external_force_resultant_N=max_force, maximum_external_reference_torque_N_m=max_torque,
        final=rows[-1], caution='Directional-Y, fixed-normal, quasi-static benchmark; not AFSD or a 2D friction cone. No independent continuum cube solution claimed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inspect(args.run)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(0 if result['status']=='PASS_COMPLETE_READBACK' else 2)
