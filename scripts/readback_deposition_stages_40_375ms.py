"""Observe the already audited 448 true frames; never rerun physical maps."""
from pathlib import Path
import hashlib
import json
import time
import traceback
import numpy as np
from deposition_stage_observer import observe_deposition_stages

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'outputs/deposition_observer_20261004/readback_40_375ms_v1'
CURRENT=ROOT/'outputs/visible_support_20261004/readback_v1'
PRIOR=ROOT/'outputs/coupled_cost_20261004/batched/readback_pre_occlusion_39_375ms'
PARENT=ROOT/'outputs/visible_support_20261004/process_v1/fine/step_016/accepted.npz'
PARENT_SHA='b03422f4fbcaaddb8151f01c397b2a6c744145f0e1f1e04a7228c879230f3dd1'
SUMMARY_SHA='54543bdc917afc3b2be431549818bf8b523815a2e36c8d48d12d6b04a95a50d8'
PRIOR_SHA='4632c3cfc7430ad00942aeea58edb227607333c5b078610ca66af1838595b600'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if OUTPUT.exists():raise ValueError('refuse to overwrite evidence')
    OUTPUT.mkdir(parents=True);start=time.perf_counter();bound={};report={};code=2
    def bind(path,expected=None):
        actual=sha(path)
        if expected is not None and actual!=expected:raise ValueError('changed source '+str(path))
        bound[path.relative_to(ROOT).as_posix()]=actual
        return path
    def read(path,expected=None):return json.loads(bind(path,expected).read_text(encoding='utf-8'))
    try:
        for path in (Path(__file__),ROOT/'scripts/deposition_stage_observer.py',ROOT/'scripts/check_deposition_stage_observer.py'):
            bind(path)
        current=read(CURRENT/'summary.json',SUMMARY_SHA);prior=read(PRIOR/'summary.json',PRIOR_SHA)
        if current['status']!='PASS_SAVED_EVIDENCE_READBACK_NOT_FRESH_PHYSICS':raise ValueError('current readback failed')
        with np.load(bind(PARENT,PARENT_SHA),allow_pickle=True) as d:
            # Read only the already hash-bound locally produced checkpoint.
            state=d['solver_data'].item();saved=d['particles'].item()['solid']['arrays']
        p=state['parameters'];events=state['bore_observation']
        with np.load(bind(ROOT/current['trajectory']['path'],current['trajectory']['sha256']),allow_pickle=False) as d:
            new={k:d[k] for k in d.files}
        with np.load(bind(PRIOR/'accepted_trajectory.npz',prior['trajectory_sha256']),allow_pickle=False) as d:
            old={k:d[k] for k in ('time_s','particle_id','position_m','plastic_strain','tool_origin_m','tool_rotation','contact_active','contact_activity_available')}
        if not np.array_equal(old['particle_id'],new['particle_id']):raise ValueError('persistent ID ordering differs')
        for key in ('time_s','position_m','plastic_strain','tool_origin_m','tool_rotation'):
            if not np.array_equal(old[key][-1],new[key][0]):raise ValueError('branch handoff differs '+key)
        args={key:np.concatenate((old[key],new[key][1:]),axis=0) for key in ('time_s','position_m','plastic_strain','tool_origin_m','tool_rotation')}
        shoulder=[old['contact_active'][:,:,0]]
        history=current['accepted_history']
        if len(history)!=16:raise ValueError('registered 16 accepted frames absent')
        for row in history:
            path=ROOT/row['checkpoint'];bind(path,row['sha256'])
            j=read(path.parent/'FULL_UNCOMMITTED_candidate.json')
            if j['contact']['features'][0]!='SHOULDER':raise ValueError('contact feature ordering differs')
            shoulder.append(np.asarray(j['contact']['active'],bool)[None,:,0])
        args.update(particle_id=new['particle_id'],component_id=new['component_id'],mass_kg=new['mass_kg'],
            eligible_ids=events['eligible_ids'],crossing_events=events['events'],dx_m=p['dx'],
            bore_width_m=p['bore_width'],tool_radius_m=p['tool_radius'],
            shoulder_contact=np.concatenate(shoulder,axis=0),
            shoulder_contact_available=np.concatenate((old['contact_activity_available'],np.ones(16,dtype=bool))))
        if len(args['time_s'])!=448 or state['step']!=447:raise ValueError('actual frame count differs')
        np.testing.assert_array_equal(args['position_m'][-1],np.column_stack([saved[k] for k in ('x','y','z')]))
        report=observe_deposition_stages(**args)
        local=(args['position_m'][-1]-args['tool_origin_m'][-1])@args['tool_rotation'][-1]
        ids=args['particle_id'];remaining=np.isin(ids,events['eligible_ids']) & ~np.isin(ids,events['crossing_ids'])
        above=np.flatnonzero(remaining & (local[:,2]>0));order=above[np.argsort(local[above,2])]
        initial=args['position_m'][0];substrate=args['component_id']==0;bar=args['component_id']==1
        report['finite_inventory_and_geometric_references']=dict(
            scope='GEOMETRY_AND_PRESCRIBED_KINEMATICS_ONLY_NOT_MATERIAL_TRANSPORT_PREDICTION',
            total_mass_kg=float(np.sum(args['mass_kg'])),bar_mass_kg=float(np.sum(args['mass_kg'][bar])),
            initial_substrate_centre_bounds_m=[initial[substrate].min(axis=0).tolist(),initial[substrate].max(axis=0).tolist()],
            initial_substrate_nominal_width_m=p['substrate_width'],tool_shoulder_diameter_m=2*p['tool_radius'],
            current_tool_origin_m=args['tool_origin_m'][-1].tolist(),
            prescribed_feed_travel_since_initial_m=p['feed']*state['time_s'],
            prescribed_traverse_travel_since_initial_m=p['traverse']*state['time_s'],
            prescribed_feed_time_per_spacing_s=p['dx']/p['feed'],
            eligible_initial_volume_kg_over_rho0_m3=report['eligible_initial_mass_kg']/p['rho0'],
            eligible_initial_equivalent_length_over_nominal_bar_area_m=report['eligible_initial_mass_kg']/p['rho0']/p['width']**2,
            nearest_uncrossed_eligible_centres_above_outlet=[dict(particle_id=int(ids[k]),height_m=float(local[k,2]),
                square_projection_margin_m=float(p['bore_width']/2-np.max(np.abs(local[k,:2])))) for k in order[:8]],
            limitations=['The finite coupon is initially 4 mm wide while the shoulder diameter is 5.4 mm.',
                'No sustained inflow is present; initial eligible mass bounds all possible new eligible crossing mass.',
                'Prescribed feed time per spacing is not the free-material transit time or an execution-time forecast.'])
        if any(sha(ROOT/path)!=expected for path,expected in bound.items()):raise ValueError('input changed during observation')
        report.update(status='PASS_PASSIVE_ACCEPTED_ARRAY_OBSERVATION',input_and_source_sha256=bound,
            accepted_checkpoint=PARENT.relative_to(ROOT).as_posix(),accepted_checkpoint_sha256=PARENT_SHA,
            actual_frames_only=True,interpolated_frames=0,added_accepted_time_s=0.)
        code=0
    except Exception:
        report.update(status='FAIL_PASSIVE_OBSERVATION',error=traceback.format_exc(),input_and_source_sha256=bound)
    report.update(exit_code=code,wall_seconds=time.perf_counter()-start)
    (OUTPUT/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report.get(k) for k in ('status','frame_count','time_end_s','eligible_initial_mass_kg','final','crossed_particle_rows','error','wall_seconds')},indent=2))
    return code


if __name__=='__main__':raise SystemExit(main())
