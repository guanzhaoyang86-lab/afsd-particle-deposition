"""Read saved checkpoints/journals only; never evaluate a physical trial."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from pysph.solver.utils import load

BASE = ROOT/'outputs/coupled_cost_20261004/batched/window_6_20ms'
TARGET = ROOT/'outputs/coupled_cost_20261004/window_error_audit.json'

def read(path):
    return json.loads(path.read_text(encoding='utf-8'))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def snapshot(branch, row):
    path = BASE/branch/row['checkpoint']
    assert digest(path) == row['sha256'], 'checkpoint hash changed'
    data = load(str(path))
    pa, state = data['arrays']['solid'], data['solver_data']
    selected = 'newton_retry' if row['fallback_reason'] else 'broyden'
    journal_path = path.parent/selected/'last_UNCOMMITTED_candidate.json'
    journal = read(journal_path)
    with np.load(journal_path.with_suffix('.npz')) as candidate:
        assert np.array_equal(np.column_stack((pa.x,pa.y,pa.z)), candidate['position_m'])
        assert np.array_equal(pa.particle_id, candidate['particle_id'])
    assert state['step'] == row['global_step']
    assert state['time_s'] == row['time_s']
    return dict(pa=pa, state=state, contact=journal['contact'],
                x=np.column_stack((pa.x,pa.y,pa.z)),
                sha256=row['sha256'], journal_sha256=digest(journal_path))

def particle_info(snap, index):
    pa, c = snap['pa'], snap['contact']
    active = np.asarray(c['active'], bool)[index]
    return dict(particle_id=int(pa.particle_id[index]), component_id=int(pa.component_id[index]),
        actuator=bool(pa.actuator[index]), active_features=[f for f,a in zip(c['features'],active) if a],
        initially_above_shoulder_eligible=bool(int(pa.particle_id[index]) in snap['state']['bore_observation']['eligible_ids']),
        sliding_features=[f for f,a in zip(c['features'],np.asarray(c['sliding'],bool)[index]) if a],
        contact_force_N=np.asarray(c['force_N'])[index].tolist(),
        position_m=snap['x'][index].tolist(), velocity_m_s=[float(pa.u[index]),float(pa.v[index]),float(pa.w[index])],
        ep=float(pa.ep[index]))

def geometry_summary(snap):
    pa, state, c = snap['pa'], snap['state'], snap['contact']
    p = state['parameters']; all_bar = np.flatnonzero(pa.component_id == 1)
    bar = np.flatnonzero(np.isin(pa.particle_id,state['bore_observation']['eligible_ids']))
    z = snap['x'][bar,2]-p['layer']; nearest = int(bar[np.argmin(z)])
    margins = np.asarray(c['gaps_m'])-p['dx']/2
    nonbore = {}
    for j in (0,1,6,7):
        order = np.argsort(margins[:,j])
        nonbore[c['features'][j]] = [dict(particle_id=int(pa.particle_id[i]),component_id=int(pa.component_id[i]),
             margin_to_penalty_shell_m=float(margins[i,j]),active=bool(np.asarray(c['active'])[i,j])) for i in order[:5]]
    return dict(nearest_eligible_bar_center_to_shoulder_m=float(z.min()),nearest_eligible_bar=particle_info(snap,nearest),
        all_bar_minimum_signed_plane_height_m=float((snap['x'][all_bar,2]-p['layer']).min()),
        eligibility_note='Only IDs initially above shoulder count as transported feed. Initial lower bar material is excluded from nearest waiting-center distance.',
        new_feature_margins=nonbore,crossing_mass_kg=float(state['crossing_mass_kg']),
        feature_active_counts=np.asarray(c['active'],bool).sum(axis=0).tolist())

coarse_rows, fine_rows = read(BASE/'coarse/accepted_history.json'), read(BASE/'fine/accepted_history.json')
coarse = [snapshot('coarse',r) for r in coarse_rows]
fine = [snapshot('fine',r) for r in fine_rows]
rows = []
for cr, cs in zip(coarse_rows, coarse):
    matches=[i for i,fr in enumerate(fine_rows) if abs(fr['time_s']-cr['time_s'])<1e-15]
    assert len(matches)==1
    fs=fine[matches[0]]
    assert np.array_equal(cs['pa'].particle_id,fs['pa'].particle_id)
    assert np.array_equal(cs['pa'].m,fs['pa'].m)
    error=cs['x']-fs['x']; norms=np.linalg.norm(error,axis=1); i=int(np.argmax(norms))
    components={str(int(comp)):float(norms[fs['pa'].component_id==comp].max()) for comp in np.unique(fs['pa'].component_id)}
    p=fs['state']['parameters']; t=fs['state']['time_s']
    relative_xy=fs['x'][i,:2]-[p['traverse']*t,0.]
    radial=relative_xy/np.linalg.norm(relative_xy); tangent=np.array([-radial[1],radial[0]])
    rows.append(dict(time_s=cr['time_s'],max_position_error_m=float(norms[i]),rms_position_error_m=float(np.sqrt(np.mean(norms**2))),
        max_axis_abs_error_m=np.abs(error).max(axis=0).tolist(),max_particle_error_xyz_m=error[i].tolist(),
        max_particle_error_radial_tangential_z_m=[float(error[i,:2]@radial),float(error[i,:2]@tangent),float(error[i,2])],
        component_max_position_error_m=components,max_particle_coarse=particle_info(cs,i),max_particle_fine=particle_info(fs,i),
        top_five=[dict(particle_id=int(fs['pa'].particle_id[j]),component_id=int(fs['pa'].component_id[j]),
            error_norm_m=float(norms[j]),error_xyz_m=error[j].tolist()) for j in np.argsort(norms)[-5:][::-1]],
        coarse_geometry=geometry_summary(cs),fine_geometry=geometry_summary(fs),
        coarse_sha256=cs['sha256'],fine_sha256=fs['sha256']))

t=np.array([r['time_s']-.006 for r in rows]); e=np.array([r['max_position_error_m'] for r in rows])
def fit(design):
    coefficients=np.linalg.lstsq(design,e,rcond=None)[0]; predicted=design@coefficients
    return dict(coefficients=coefficients.tolist(),rmse_m=float(np.sqrt(np.mean((predicted-e)**2))),
        r_squared=float(1-np.sum((predicted-e)**2)/np.sum((e-e.mean())**2)))
growth=dict(increments_per_1ms_m=np.diff(np.r_[0.,e]).tolist(),monotone=bool(np.all(np.diff(e)>=0.)),
    through_origin_linear=fit(t[:,None]),through_origin_linear_plus_quadratic=fit(np.column_stack((t,t*t))),
    caveat='Two step sizes cannot establish an observed convergence order. Growth fit is descriptive only.',
    conditional_next_halving_error_if_first_order_m=float(e[-1]/2),
    fixed_position_gate_m=5e-6)

# Pure extrapolation of SAVED per-particle distances, not new tool/contact calls.
window=[i for i,r in enumerate(fine_rows) if r['time_s']>=.015-1e-15]
times=np.array([fine_rows[i]['time_s'] for i in window]); last=fine[-1]; p=last['state']['parameters']
gaps=np.stack([np.asarray(fine[i]['contact']['gaps_m']) for i in window]); centered=times-times.mean()
slopes=np.sum(centered[:,None,None]*gaps,axis=0)/np.sum(centered**2)
margins=gaps[-1]-p['dx']/2
predictions=[]
for j in (0,1,6,7):
    for i in range(len(last['pa'].x)):
        if margins[i,j]>1e-12 and slopes[i,j]<-1e-8:
            extra=float(-margins[i,j]/slopes[i,j])
            predictions.append(dict(feature=last['contact']['features'][j],particle_id=int(last['pa'].particle_id[i]),
                component_id=int(last['pa'].component_id[i]),margin_m=float(margins[i,j]),
                saved_distance_linear_slope_m_s=float(slopes[i,j]),extrapolated_additional_time_s=extra,
                extrapolated_total_time_s=float(times[-1]+extra)))
predictions.sort(key=lambda r:r['extrapolated_additional_time_s'])
bar=np.flatnonzero(last['pa'].component_id==1); positions=np.stack([fine[i]['x'] for i in window])
zslopes=np.sum(centered[:,None]*positions[:,:,2],axis=0)/np.sum(centered**2)
cross_predictions=[]
for i in bar:
    height=last['x'][i,2]-p['layer']
    if height>0 and zslopes[i]<-1e-8:
        cross_predictions.append(dict(particle_id=int(last['pa'].particle_id[i]),height_m=float(height),
            saved_z_slope_m_s=float(zslopes[i]),extrapolated_total_time_s=float(times[-1]-height/zslopes[i])))
cross_predictions.sort(key=lambda r:r['extrapolated_total_time_s'])
report=dict(scope='READ_ONLY_SAVED_CHECKPOINT_AND_JOURNAL_ANALYSIS_NO_PHYSICAL_EVALUATION',
    controller_source_sha256=digest(Path(__file__)),coarse_status=read(BASE/'coarse/manifest.json')['status'],
    fine_status=read(BASE/'fine/manifest.json')['status'],
    endpoint_comparison=read(BASE/'fine/temporal_comparison.json'),equal_time_rows=rows,error_growth=growth,
    final_fine_geometry=geometry_summary(last),
    saved_geometry_event_extrapolation=dict(fit_interval_s=[float(times[0]),float(times[-1])],
        earliest_nonbore_shell_candidates=predictions[:12],earliest_bar_center_plane_candidates=cross_predictions[:8],
        caveat='Unvalidated linear extrapolation of saved gaps only; no claim of first event time, path safety, crossing or transport. Distances and rates can change, and substep shell events remain unresolved by sparse samples.'))
refinement_comparisons=[]
previous=last; previous_error=float(e[-1])
for level,count,step in (('quarter',56,.00025),('eighth',112,.000125)):
    branch='refinement/'+level
    history_path=BASE/branch/'accepted_history.json'
    if not history_path.exists():break
    saved_rows=read(history_path)
    if len(saved_rows)!=count or abs(saved_rows[-1]['time_s']-.020)>1e-15:break
    snap=snapshot(branch,saved_rows[-1]);diff=previous['x']-snap['x'];norms=np.linalg.norm(diff,axis=1);i=int(norms.argmax())
    current_error=float(norms[i])
    refinement_comparisons.append(dict(level=level,dt_s=step,position_error_against_previous_m=current_error,
        max_particle=particle_info(snap,i),coarser_minus_finer_xyz_m=diff[i].tolist(),
        previous_error_ratio=float(previous_error/current_error),observed_pair_order=np.log2(previous_error/current_error),
        geometry=geometry_summary(snap),sha256=snap['sha256'],
        interpretation='Ratio of successive maximum-position differences; same PID required before interpreting as particlewise order. No exact solution comparison.'))
    previous=snap;previous_error=current_error
report['later_completed_refinement_endpoints']=refinement_comparisons
if TARGET.exists():raise ValueError('refuse overwrite')
TARGET.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps(dict(path=str(TARGET),max_rows=[dict(t_ms=r['time_s']*1e3,error_um=r['max_position_error_m']*1e6,pid=r['max_particle_fine']['particle_id'],component=r['max_particle_fine']['component_id'],active=r['max_particle_fine']['active_features']) for r in rows],growth=growth,final_geometry=report['final_fine_geometry'],event_estimates=predictions[:4],crossing_estimates=cross_predictions[:2]),indent=2))
