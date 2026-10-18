"""Read-only inertial-force scale on the accepted 0--30 ms QS trajectory.

The saved velocity is verified against each actual displacement/dt. It is a
step-average kinematic rate, not a dynamically solved instantaneous velocity.
At shared node t[n], use 2*(vbar[n+1]-vbar[n])/(dt[n+1]+dt[n]). This three-point
second derivative uses the real unequal intervals and is exact for quadratic
position histories. Compare m*a with separately saved internal/contact force
constituents at that same actual node, never with their near-zero balance.

No trajectory interpolation, trial construction, physical evaluation, contact
geometry, NNPS, or material update is performed. There is no validity threshold.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import readback_discrete_process as reader

ROOT = Path(__file__).resolve().parents[1]
READBACK = ROOT/'outputs/coupled_cost_20261004/batched/readback_30ms_refined'
DEFAULT_OUTPUT = ROOT/'outputs/coupled_cost_20261004/inertia_omission_30ms'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def acceleration(position, times):
    dt = np.diff(times)
    require(np.all(dt > 0.), 'strictly increasing actual sample times required')
    average_velocity = np.diff(position, axis=0)/dt[:, None, None]
    a = np.diff(average_velocity, axis=0)/(.5*(dt[:-1]+dt[1:]))[:, None, None]
    return average_velocity, a


def ratio(numerator, denominator):
    return None if denominator == 0. else float(numerator/denominator)


def statistics(values):
    values = np.asarray([value for value in values if value is not None], dtype=float)
    if not len(values):
        return dict(count=0, minimum=None, median=None, p95=None, maximum=None)
    require(np.isfinite(values).all(), 'nonfinite diagnostic quantity')
    return dict(count=len(values), minimum=float(values.min()), median=float(np.median(values)),
                p95=float(np.quantile(values, .95)), maximum=float(values.max()))


def run(output):
    begin = time.perf_counter()
    output = Path(output).resolve()
    json_path, arrays_path = output.with_suffix('.json'), output.with_suffix('.npz')
    require(not json_path.exists() and not arrays_path.exists(), 'refuse diagnostic evidence overwrite')
    inputs = {}

    def bind(path, expected=None):
        digest = reader.sha(path)
        require(expected is None or expected == digest, 'bound input SHA mismatch: '+str(path))
        inputs[path.relative_to(ROOT).as_posix()] = digest
        return path

    summary = reader.read(bind(READBACK/'summary.json'))
    require(summary['status'].startswith('PASS_') and summary['accepted_steps']==281 and summary['samples']==282
            and abs(summary['accepted_time_s']-.030)<1e-15, 'complete accepted 30 ms readback required')
    bind(READBACK/'accepted_trajectory.npz', summary['trajectory_sha256'])
    bound_readback_inputs = summary['source_and_input_sha256']
    with np.load(READBACK/'accepted_trajectory.npz', allow_pickle=False) as saved:
        require(bool(saved['actual_checkpoint_samples_only']) and not bool(saved['interpolation_performed']),
                'only actual checkpoint samples supported')
        keys = ('time_s','particle_id','position_m','velocity_m_s','deviatoric_stress_Pa','temperature_K',
                'density_kg_m3','plastic_strain','contact_history_m','particle_property__m',
                'particle_property__component_id')
        fields = {key:saved[key].copy() for key in keys}
    times, ids, position = fields['time_s'], fields['particle_id'], fields['position_m']
    require(len(times)==282 and np.array_equal(ids,np.sort(np.unique(ids))), 'actual sample shape/order differs')
    mass_frames = fields['particle_property__m'][...,0]
    require(np.array_equal(mass_frames,np.broadcast_to(mass_frames[0],mass_frames.shape)) and np.all(mass_frames>0),
            'constant positive particle masses required')
    mass = mass_frames[0]
    initial_path = ROOT/'outputs/discrete_contact_process_20261004/v3/quarters/independent_initial.npz'
    initial, initial_state = reader.restore(bind(initial_path,bound_readback_inputs[initial_path.relative_to(ROOT).as_posix()]))
    order_initial = np.argsort(initial.particle_id)
    require(np.array_equal(initial.particle_id[order_initial],ids), 'initial free-mask ID mapping differs')
    free = np.asarray(initial_state['free'],bool)[order_initial]
    require(free.shape==(len(ids),3), 'expected actual physical free-DOF mask')
    require(np.array_equal(fields['velocity_m_s'][0],np.zeros((len(ids),3))), 'independent initial velocity is not zero')
    average_velocity, a = acceleration(position,times)
    dt = np.diff(times)
    velocity_error = np.linalg.norm(fields['velocity_m_s'][1:]-average_velocity,axis=2)
    require(float(velocity_error.max())<1e-10, 'saved velocity is not the displacement/dt rate')
    # Pure arithmetic sanity check of the unequal-dt formula, no model test.
    fixture_times = np.array([0.,.13,.4,.71,1.2])
    fixture_position = (1.+3.*fixture_times+2.*fixture_times**2)[:,None,None]*np.ones((1,1,3))
    require(np.allclose(acceleration(fixture_position,fixture_times)[1],4.,rtol=0.,atol=2e-13),
            'unequal-time quadratic acceleration check failed')
    force_internal, force_contact, journal_rows = [], [], []
    checks = (('position_m','position_m',1e-15),('velocity_m_s','velocity_m_s',1e-12),
              ('stress_Pa','deviatoric_stress_Pa',1e-6),('temperature_K','temperature_K',1e-9),
              ('density_kg_m3','density_kg_m3',1e-9),('plastic_strain','plastic_strain',1e-14),
              ('contact_history','contact_history_m',1e-15))
    for index,row in enumerate(summary['rows'],1):
        checkpoint = bind(ROOT/row['checkpoint'],row['sha256'])
        require(row['step']==index and row['time_s']==times[index] and abs(row['dt_s']-dt[index-1])<1e-15,
                'accepted row time/step differs from trajectory')
        matches = []
        for path in sorted(checkpoint.parent.rglob('last_UNCOMMITTED_candidate.npz')):
            name = path.relative_to(ROOT).as_posix()
            if name not in bound_readback_inputs:
                continue
            with np.load(path,allow_pickle=False) as saved:
                order = np.argsort(saved['particle_id'])
                if not np.array_equal(saved['particle_id'][order],ids):
                    continue
                if not all(np.allclose(saved[key][order],fields[target][index],rtol=0.,atol=tol)
                           for key,target,tol in checks):
                    continue
                fi,fc = saved['internal_force_N'][order].copy(),saved['contact_force_N'][order].copy()
                require(np.isfinite(fi).all() and np.isfinite(fc).all(), 'nonfinite saved force constituent')
                require(np.allclose(saved['force_N'][order],fi+fc,rtol=1e-13,atol=1e-13),
                        'force journal sum does not match its separately saved constituents')
            bind(path,bound_readback_inputs[name]);matches.append((path,fi,fc))
        require(len(matches)==1, 'expected one readback-bound matching successful force journal: '+str(checkpoint))
        path,fi,fc=matches[0];force_internal.append(fi);force_contact.append(fc)
        journal_rows.append(dict(step=index,time_s=times[index],stage=row['stage'],
            checkpoint=row['checkpoint'],checkpoint_sha256=row['sha256'],
            force_journal=path.relative_to(ROOT).as_posix(),force_journal_sha256=reader.sha(path)))
    fi_all,fc_all=np.array(force_internal),np.array(force_contact)
    # a[j] is associated with accepted shared node j+1, not the terminal node.
    fi,fc=fi_all[:-1],fc_all[:-1]
    inertia=mass[None,:,None]*a
    components=fields['particle_property__component_id'][0,:,0]
    groups=dict(all_dofs=np.ones_like(free,bool),physical_free_dofs=free,
        bar_free_dofs=free&(components==1)[:,None],substrate_free_dofs=free&(components==0)[:,None])
    rows=[]
    for n in range(len(a)):
        forces={}
        for name,mask in groups.items():
            inertial_norm=float(np.linalg.norm(inertia[n][mask]))
            internal_norm=float(np.linalg.norm(fi[n][mask]));contact_norm=float(np.linalg.norm(fc[n][mask]))
            scale=float(np.hypot(internal_norm,contact_norm))
            forces[name]=dict(inertial_proxy_l2_N=inertial_norm,internal_force_l2_N=internal_norm,
                contact_force_l2_N=contact_norm,constituent_force_scale_l2_N=scale,
                inertial_to_constituent_ratio=ratio(inertial_norm,scale))
        particle_inertia=np.linalg.norm(inertia[n]*free,axis=1)
        particle_scale=np.sqrt(np.sum((fi[n]*free)**2+(fc[n]*free)**2,axis=1))
        # Descriptive subset only, not an acceptance cutoff. Keep raw values too.
        support_floor=.01*float(particle_scale.max())
        supported=(particle_scale>0.)&(particle_scale>=support_floor)
        supported_ratio=particle_inertia[supported]/particle_scale[supported]
        worst=int(np.argmax(particle_inertia))
        rows.append(dict(step=n+1,time_s=float(times[n+1]),stage=summary['rows'][n]['stage'],
            left_dt_s=float(dt[n]),right_dt_s=float(dt[n+1]),
            dt_changed=not bool(np.isclose(dt[n],dt[n+1],rtol=1e-12,atol=1e-15)),
            acceleration_max_m_s2=float(np.max(np.linalg.norm(a[n],axis=1))),
            groups=forces,maximum_particle_free_inertial_proxy_N=float(particle_inertia.max()),
            particle_at_max_free_inertial_proxy=int(ids[worst]),
            constituent_scale_at_max_inertial_particle_N=float(particle_scale[worst]),
            supported_pointwise_scale_floor_N=support_floor,supported_particle_count=int(supported.sum()),
            unsupported_or_small_scale_particle_count=int((~supported).sum()),
            supported_pointwise_ratio_statistics=statistics(supported_ratio)))
    stage_stats={}
    for stage in sorted({row['stage'] for row in rows}):
        subset=[row for row in rows if row['stage']==stage]
        stage_stats[stage]=dict(nodes=len(subset),start_s=subset[0]['time_s'],end_s=subset[-1]['time_s'],
            free_ratio_statistics=statistics([r['groups']['physical_free_dofs']['inertial_to_constituent_ratio'] for r in subset]))
    aggregate={name:statistics([row['groups'][name]['inertial_to_constituent_ratio'] for row in rows]) for name in groups}
    free_ratios=[row['groups']['physical_free_dofs']['inertial_to_constituent_ratio'] for row in rows]
    worst_index=max(range(len(rows)),key=lambda n:free_ratios[n] if free_ratios[n] is not None else -1.)
    initial_jump=average_velocity[0]-fields['velocity_m_s'][0]
    startup=dict(initial_velocity_is_initialized_not_dynamically_solved=True,
        first_interval_s=float(dt[0]),first_average_velocity_time_s=float(.5*(times[0]+times[1])),
        initial_to_first_average_velocity_jump_max_m_s=float(np.max(np.linalg.norm(initial_jump,axis=1))),
        mass_times_velocity_jump_l2_kg_m_s=float(np.linalg.norm(mass[:,None]*initial_jump)),
        actual_startup_acceleration=None,
        reason='No resolved startup ramp or instantaneous endpoint velocities; acceleration of an abrupt prescribed start is not identified by these samples.')
    np.savez_compressed(arrays_path,particle_id=ids,mass_kg=mass,physical_free_dofs=free,
        step_average_velocity_time_s=.5*(times[:-1]+times[1:]),step_average_velocity_m_s=average_velocity,
        acceleration_node_time_s=times[1:-1],left_dt_s=dt[:-1],right_dt_s=dt[1:],
        derived_acceleration_m_s2=a,inertial_force_proxy_N=inertia,
        saved_internal_force_N=fi,saved_contact_force_N=fc,
        step_average_velocity_kinetic_scale_J=.5*np.sum(mass[None,:]*np.sum(average_velocity**2,axis=2),axis=1),
        physical_solver_evaluations=np.array(0),diagnostic_not_dynamic_validation=np.array(True))
    source_paths=(Path(__file__).resolve(),Path(reader.__file__).resolve(),
        ROOT/'src/afsd_backend/pysph/quasistatic_process.py',ROOT/'src/afsd_backend/pysph/discrete_process.py')
    report=dict(status='COMPUTED_READ_ONLY_INERTIAL_SCALE_NOT_DYNAMIC_OR_INERTIA_VALIDATION',
        physical_solver_evaluations=0,nnps_evaluations=0,material_updates=0,
        accepted_endpoint_time_s=float(times[-1]),accepted_steps=281,actual_samples=282,
        acceleration_node_count=len(a),last_acceleration_node_s=float(times[-2]),
        terminal_acceleration_available=False,criterion_or_timestep_changed=False,
        velocity_semantics='Saved QS displacement/dt, verified against actual accepted positions; interval-average kinematic rate.',
        saved_velocity_vs_displacement_rate_max_error_m_s=float(velocity_error.max()),
        acceleration_definition='a[n]=2*((x[n+1]-x[n])/dt_right-(x[n]-x[n-1])/dt_left)/(dt_left+dt_right).',
        force_denominator='sqrt(sum internal_force^2 + sum contact_force^2), on the declared DOF subset; no net residual in denominator.',
        pointwise_diagnostic_subset='Particles with constituent force scale >=1% of that node maximum; omitted small-scale particles counted, not asserted safe.',
        arithmetic_check='PASS_UNEQUAL_DT_QUADRATIC_POSITION_SECOND_DERIVATIVE',startup=startup,
        aggregate_ratios=aggregate,stage_ratios=stage_stats,worst_free_ratio_node=rows[worst_index],
        dt_transition_nodes=[row for row in rows if row['dt_changed']],
        constant_dt_free_ratio_statistics=statistics([row['groups']['physical_free_dofs']['inertial_to_constituent_ratio'] for row in rows if not row['dt_changed']]),
        rows=rows,force_journal_rows=journal_rows,
        limitations=['No comparison with a resolved inertial solve; this is a consistency scale of the existing QS path.',
            'Differencing QS motion can contain equilibrium/step-size noise and does not recover an unrecorded startup ramp.',
            'Constituent nodal internal forces may themselves cancel pairwise stresses; no per-pair stress-force decomposition is available here.',
            'Global small ratios alone cannot establish dynamic stability, local accuracy, kernel consistency, or experimental agreement.',
            'No acceleration assigned to final 30ms endpoint without a subsequent accepted sample; initial startup acceleration left unknown.'],
        input_sha256=inputs,arrays_sha256=reader.sha(arrays_path),
        source_sha256={path.relative_to(ROOT).as_posix():reader.sha(path) for path in source_paths},
        wall_seconds=time.perf_counter()-begin)
    json_path.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(status=report['status'],output=str(json_path),nodes=len(a),
        velocity_error_m_s=report['saved_velocity_vs_displacement_rate_max_error_m_s'],
        free_ratio_statistics=aggregate['physical_free_dofs'],worst_free_time_s=rows[worst_index]['time_s'],
        startup=startup,wall_seconds=report['wall_seconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    run(parser.parse_args().output)
