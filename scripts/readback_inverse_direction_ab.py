"""Independent saved-array/algebra audit; no solver, tangent, NNPS or geometry calls."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.interval_acceptance import (
    validate_interval_energy_state,validate_bore_crossing_state,observe_square_bore_crossings)

PARENT_SHA='12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b'
CLASSIFICATION='SAME_PARENT_VERIFIED_INVERSE_NUMERICAL_DIAGNOSTIC_NOT_ACCEPTED_V1'
PROCESS_CLASS='EXPANDED_HOT_QS_VISIBLE_SUPPORT_APPROXIMATE_V1'
DRIVER_PASS='PASS_SAME_PARENT_INVERSE_NUMERICAL_SCREEN_NOT_PROCESS_ACCEPTANCE'
DT=7.8125e-6
LIMITS=dict(position_max_m=1e-9,temperature_max_K=1e-4,stress_pointwise_relative=1e-3,
    stress_floor_Pa=1e6,plastic_absolute=1e-7,density_relative=1e-6,contact_history_max_m=1e-9,
    energy_absolute_J=1e-12,energy_relative=1e-3)
IMMUTABLE=('particle_id','x0','y0','z0','m','h','material_id','component_id','tracer_type','wholebar_tag','actuator')


def require(value,message):
    if not bool(value):raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def project(path):
    p=Path(path);p=(ROOT/p).resolve() if not p.is_absolute() else p.resolve()
    require(p.is_relative_to(ROOT),'input path escaped project');return p


def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    return value


def write(path,value):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(clean(value),f,indent=2,allow_nan=False);f.write('\n')


def arrays(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}


def checkpoint(path):
    # Only locally hash-bound PySPH v2 archives; never an untrusted network NPZ.
    with np.load(path,allow_pickle=True) as z:
        require(int(z['version'])==2,'unsupported serialized checkpoint version')
        particles=z['particles'].item();state=z['solver_data'].item()
        require(set(particles)=={'solid'},'unexpected particle array set')
        pa={k:np.asarray(v).copy() for k,v in particles['solid']['arrays'].items()}
    return pa,state


def pose(p,t):
    c,s=np.cos(p['omega']*t),np.sin(p['omega']*t)
    return dict(origin=np.array([p['traverse']*t,0.,p['layer']]),rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]))


def values_difference(a,b):
    require(set(a)==set(b),'saved particle property set differs')
    rows={}
    for key in sorted(a):
        x,y=np.asarray(a[key]),np.asarray(b[key])
        require(x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all(),'invalid saved field '+key)
        d=x.astype(float)-y.astype(float)
        rows[key]=dict(shape=list(x.shape),equal=bool(np.array_equal(x,y)),
            maximum_absolute_difference=float(np.max(np.abs(d))) if d.size else 0.,l2_difference=float(np.linalg.norm(d.ravel())))
    return rows


class Reader:
    def __init__(self,directory):self.directory=directory;self.evidence={}

    def bind(self,path,digest=None):
        path=project(path);actual=sha(path)
        require(digest is None or digest==actual,'source/artifact hash drift: '+rel(path))
        key=rel(path);require(key not in self.evidence or self.evidence[key]==actual,'input changed while reading: '+key)
        self.evidence[key]=actual;return path

    def mapping(self,mapping):
        for name,digest in mapping.items():self.bind(name,digest)

    def initial(self,parent,state):
        folder=self.directory/'initial_direction'
        z=arrays(self.bind(folder/'same_initial_linear_system_NOT_RESTARTABLE.npz'))
        journal=arrays(self.bind(folder/'INITIAL_UNCOMMITTED_NOT_EQUILIBRIUM.npz'))
        metrics=read(self.bind(folder/'comparison.json'))
        A,H,rhs=z['A'],z['H'],z['rhs'];free=z['free'].astype(bool)
        require(A.shape==H.shape==(1961,1961) and rhs.shape==(1961,) and free.shape==(918,3),
            'wrong saved common linear system shape')
        require(np.array_equal(free,state['free']) and np.array_equal(z['particle_id'],parent['particle_id']) and
            not bool(z['accepted']) and not bool(z['restart_allowed']),'initial basis or classification differs')
        guess=np.column_stack([parent[k] for k in ('u','v','w')])*DT
        prescribed=np.zeros_like(guess);prescribed[parent['actuator'].astype(bool),2]=-state['parameters']['feed']*DT
        guess[~free]=prescribed[~free]
        require(np.array_equal(z['guess_m'],guess) and np.array_equal(journal['displacement_m'],guess) and
            np.array_equal(journal['particle_id'],parent['particle_id']) and np.array_equal(rhs,-journal['force_N'][free]),
            'common system did not use actual parent guess / row ordering / force rhs')
        condition=float(np.linalg.norm(A,np.inf)*np.linalg.norm(H,np.inf))
        limit=1e-9*max(1.,float(np.max(np.abs(rhs))))
        errors={key:float(np.max(np.abs(A@z[key]-rhs))) for key in ('original_direction_m','inverse_direction_m')}
        require(all(np.isfinite(v) and v<=limit for v in errors.values()) and condition<=1e12,
            'saved initial directions fail the original matrix/backward gate')
        require(abs(condition-metrics['condition_infinity_estimate'])<=1e-10*max(condition,1.),'initial condition log differs')
        # Both directions are verified against the one saved A/rhs, never two near-identical systems.
        delta=float(np.max(np.abs(z['original_direction_m']-z['inverse_direction_m'])))
        require(np.isclose(delta,metrics['direction_max_absolute_difference_m'],rtol=1e-10,atol=1e-25),
            'initial direction difference log differs')
        return dict(status='PASS_ACTUAL_SAVED_COMMON_MATRIX_RHS_AND_BOTH_BACKWARD_ERRORS',
            matrix_file_sha256=sha(folder/'same_initial_linear_system_NOT_RESTARTABLE.npz'),
            dimension=1961,condition_infinity_estimate=condition,linear_backward_limit_N=limit,
            recomputed_backward_errors_N=errors,direction_max_difference_m=delta,
            reported_timing=metrics,initial_tangent_reassembled=False)

    def iteration_logs(self,folder,route):
        checks=[];directions=[];refreshes=0;methods={};map_times=[]
        for attempt in ('broyden','newton_retry'):
            path=folder/attempt/'iterations.json'
            if not path.exists():continue
            rows=read(self.bind(path));require(rows and len(rows)<=25,'iteration budget / log count differs')
            require(max(r['evaluations'] for r in rows)<=100,'recorded evaluate budget exceeded')
            if attempt=='broyden':
                for row in rows:
                    require(np.isclose(row['actual_free_force_residual_N'],row['force_residual_N'],rtol=0.,atol=0.),
                        'actual and logged nonlinear force residual differ')
                    for build in row.get('jacobian_builds',[]):
                        require(build['condition_infinity']<=1e12,'recorded tangent condition exceeded')
                        refreshes+=int(build['build_number']>1)
                    for lin in row.get('linear',[]):
                        # ||rhs||inf is the saved actual free-force maximum at this iterate.
                        limit=1e-9*max(1.,row['actual_free_force_residual_N'])
                        require(lin['linear_backward_error_N']<=limit and lin['fixed_axis_increment_zero'],
                            'recorded per-direction backward gate failed')
                        item=dict(iteration=row['iteration'],linear_backward_error_N=lin['linear_backward_error_N'],
                            limit_N=limit,actual_A_d_not_resaved=True)
                        if route=='inverse':
                            v=lin['verified_inverse_direction'];method=v['method'];methods[method]=methods.get(method,0)+1
                            require(v['linear_backward_error_N']==lin['linear_backward_error_N'] and
                                v['linear_backward_limit_N']==limit and v['condition_infinity_estimate']<=1e12 and
                                v['matrix_size']==1961,'inverse direction diagnostics differ')
                            require(v['correction_beneficial'] is False or
                                (v['correction_attempted'] and v['correction_backward_error_N']<v['inverse_initial_backward_error_N']),
                                'non-improving correction was kept')
                            require(v['dense_fallback_used']==(method=='VERIFIED_ORIGINAL_DENSE_FALLBACK'),
                                'dense fallback route flag differs')
                            item['verified_inverse']=v
                        directions.append(item)
                    secant=row.get('secant',{})
                    if 'condition_infinity_estimate' in secant:require(secant['condition_infinity_estimate']<=1e12,'secant condition exceeded')
                require(refreshes<=2,'recorded refresh budget exceeded')
            map_times.extend(r['wall_seconds'] for r in rows)
            checks.append(dict(attempt=attempt,stored_iterates=len(rows),last_evaluations=rows[-1]['evaluations'],
                terminal_force_residual_N=rows[-1]['force_residual_N'],terminal_force_limit_N=rows[-1]['force_limit_N']))
        require(checks,'missing root attempt logs')
        return dict(attempts=checks,directions=directions,inverse_route_counts=methods,
            saved_iterate_map_seconds_sum=sum(map_times),all_evaluation_times_available=False,
            limitation='Only the initial common A/H/rhs/directions were saved; later A@d cannot be independently recomputed.')

    def candidate(self,route,row,parent,state,protocol):
        folder=self.directory/route;path=folder/'NOT_RESTARTABLE_candidate.npz'
        require(row['route']==route and row['parent_sha256']==PARENT_SHA and row['dt_s']==DT and row['global_step']==15 and
            self.directory/row['checkpoint']==path,'candidate row source/time/direct parent mismatch')
        for name,digest in row['artifacts'].items():
            p=project(self.directory/name);require(p.is_relative_to(folder),'candidate artifact outside its route');self.bind(p,digest)
        self.bind(path,row['sha256']);pa,cs=checkpoint(path)
        require(cs['accepted'] is False and cs['restart_allowed'] is False and cs['diagnostic_only'] is True and
            cs['diagnostic_classification']==CLASSIFICATION and cs['diagnostic_parent_sha256']==PARENT_SHA and
            cs['continuation_parent_sha256']==PARENT_SHA and cs['window_origin_sha256']==PARENT_SHA and
            cs['classification']==PROCESS_CLASS and cs['branch_origin_sha256']==state['branch_origin_sha256'] and
            cs['parameters']==state['parameters'] and cs['source_sha256']==state['source_sha256'] and
            cs['implementation_source_sha256']==protocol['source_sha256'] and cs['controller_source_sha256']==protocol['source_sha256'],
            'unaccepted independent candidate source / origin / branch flags differ')
        require(cs['time_s']==row['time_s']==protocol['candidate_time_s'] and cs['step']==15 and
            cs['interval_energy']['steps']==15 and cs['interval_energy']['records'][:-1]==state['interval_energy']['records'],
            'candidate endpoint or complete 14-record inherited prefix differs')
        validate_interval_energy_state(cs['interval_energy']);validate_bore_crossing_state(cs['bore_observation'])
        for key in IMMUTABLE:require(np.array_equal(pa[key],parent[key]),'immutable candidate field changed '+key)
        require(np.array_equal(cs['free'],state['free']),'candidate free axes changed')
        full=arrays(self.bind(folder/'FULL_UNCOMMITTED_candidate.npz'));detail=read(self.bind(folder/'FULL_UNCOMMITTED_candidate.json'))
        require(not bool(full['accepted']) and not bool(full['restart_allowed']) and
            detail['accepted'] is False and detail['restart_allowed'] is False,'full map diagnostic flags differ')
        for key,actual in dict(particle_id=pa['particle_id'],position_m=np.column_stack([pa[k] for k in ('x','y','z')]),
                velocity_m_s=np.column_stack([pa[k] for k in ('u','v','w')]),density_kg_m3=pa['rho'],
                stress_Pa=pa['S'].reshape(-1,3,3),temperature_K=pa['temperature'],plastic_strain=pa['ep'],
                contact_history=cs['contact_history']).items():
            require(np.array_equal(full[key],actual),'checkpoint and fresh full-map output differ '+key)
        require(np.array_equal(full['velocity_m_s'],full['displacement_m']/DT),'candidate velocity no longer actual interval displacement')
        force=float(np.max(np.abs(full['force_N'][state['free']])))
        limit=float(detail['metrics']['force_limit_N'])
        require(force<=limit and force==detail['metrics']['force_residual_N'] and detail['metrics']['scaled_residual']<=1.,
            'fresh actual saved force fails original equilibrium gate')
        energy=read(self.bind(folder/'energy.json'))
        require(energy['interval']==cs['interval_energy'] and energy['mechanical']==detail['energy']==clean(cs['energy']) and
            energy['thermal']==detail['thermal'],'saved full-map energy/thermal and checkpoint differ')
        require(energy['interval']['records'][-1]=={k:energy['mechanical'][k] for k in ('global_raw_residual_J','absolute_work_scale_J')},
            'latest full-map energy not the appended prefix record')
        geometry=read(self.bind(folder/'geometry_audit.json'))
        require(not geometry['transient_contacts_missed_by_endpoints'] and
            set(geometry['new_non_bore_features'])<={'SHOULDER','TOP'},'saved geometry outside qualified domain')
        require(read(self.bind(folder/'continuous_shell_audit.json'))['status']=='CERTIFIED_NO_HIDDEN_SHELL_CONTACT',
            'saved shell certificate incomplete')
        path_certificate=read(self.bind(folder/'visible_path_audit.json'))
        require(path_certificate['status']=='CERTIFIED_VISIBLE_SUPPORT_PATH' and path_certificate['candidate_sha256']==sha(path) and
            path_certificate['graph_sha256']==sha(folder/'visible_graph.npz'),'path certificate not for this candidate graph')
        require(read(self.bind(folder/'scalar_batch_equivalence.json'))['status']=='PASS_SAME_FILTERED_COMPLETE_MAP',
            'saved same-map complete qualification absent')
        graph=arrays(self.bind(folder/'visible_graph.npz'))
        require(np.array_equal(graph['particle_id'],pa['particle_id']) and
            np.array_equal(graph['position1'],full['position_m']),'saved actual graph row basis/position differs')
        previous=np.column_stack([parent[k] for k in ('x','y','z')])
        observed=observe_square_bore_crossings(state['bore_observation'],previous,full['position_m'],
            pa['particle_id'],pa['component_id'],pa['m'],state['time_s'],cs['time_s'],lambda t:pose(state['parameters'],t))
        require(observed['next_state']==cs['bore_observation'],'saved strict crossing observer replay differs')
        return dict(pa=pa,state=cs,full=full,detail=detail,graph=graph,row=row,energy=energy,
            logs=self.iteration_logs(folder,route),force_residual_N=force,force_limit_N=limit)

    def compare(self,a,b,declared):
        aa,bb=a['full'],b['full'];n=len(aa['particle_id'])
        sdiff=np.linalg.norm((aa['stress_Pa']-bb['stress_Pa']).reshape(n,-1),axis=1)
        srel=sdiff/np.maximum(np.linalg.norm(aa['stress_Pa'].reshape(n,-1),axis=1),1e6)
        metrics=dict(position_max_m=float(np.linalg.norm(aa['position_m']-bb['position_m'],axis=1).max()),
            temperature_max_K=float(np.abs(aa['temperature_K']-bb['temperature_K']).max()),
            stress_pointwise_relative=float(srel.max()),plastic_absolute=float(np.abs(aa['plastic_strain']-bb['plastic_strain']).max()),
            density_relative=float(np.max(np.abs(aa['density_kg_m3']-bb['density_kg_m3'])/np.maximum(np.abs(aa['density_kg_m3']),1.))),
            contact_history_max_m=float(np.linalg.norm((aa['contact_history']-bb['contact_history']).reshape(-1,3),axis=1).max()))
        for key,value in metrics.items():require(np.isclose(value,declared['metrics'][key],rtol=1e-12,atol=1e-24),
            'reported numerical error differs from saved arrays '+key)
        branches={}
        for key,value in a['detail']['contact'].items():
            if isinstance(value,list):
                arr=np.asarray(value)
                if arr.dtype.kind in 'biu':branches[key]=np.array_equal(arr,np.asarray(b['detail']['contact'][key]))
        require(branches==declared['contact_integer_boolean_branch_match'],'saved branch comparison differs')
        graph_same=all(np.array_equal(a['graph'][k],b['graph'][k]) for k in ('particle_id','i','j','raw_i','raw_j'))
        graph_same=graph_same and a['row']['graph_fingerprint']==b['row']['graph_fingerprint']
        energy=[]
        for key,value in a['energy']['mechanical'].items():
            if key.endswith('_J') and np.ndim(value)==0:
                other=b['energy']['mechanical'][key];limit=1e-12+.001*max(abs(value),abs(other))
                energy.append(dict(field=key,absolute_difference_J=abs(value-other),limit_J=limit,passed=abs(value-other)<=limit))
        crosses=a['state']['bore_observation']['crossing_ids']==b['state']['bore_observation']['crossing_ids']
        passed=all(metrics[k]<=LIMITS[k] for k in metrics) and all(branches.values()) and graph_same and crosses and all(e['passed'] for e in energy)
        require(passed==declared['passed'],'saved comparison status does not follow the preregistered limits')
        return dict(passed=passed,metrics=metrics,limits=LIMITS,graph_topology_match=graph_same,crossing_ids_match=crosses,
            branch_flags=branches,energy_comparison=energy,maximum_stress_error_particle_id=int(aa['particle_id'][np.argmax(srel)]),
            particle_property_differences=values_difference(a['pa'],b['pa']),
            full_journal_differences=values_difference(aa,bb))

    def run(self):
        d=self.directory;manifest=read(self.bind(d/'manifest.json'));protocol=read(self.bind(d/'protocol.json'))
        receipt=read(self.bind(d/'supervision.json'))
        require(manifest['status']!='RUNNING_NUMERICAL_DIAGNOSTIC' and receipt['manifest_sha256']==sha(d/'manifest.json') and
            receipt['log_sha256']==sha(d/'run.log') and receipt['source_inputs_still_match'] and
            receipt['wall_seconds']<=605.,'completed externally supervised A/B required')
        self.bind(d/'run.log');self.mapping(protocol['source_sha256']);self.mapping(protocol['reference_sha256']);self.mapping(manifest['artifacts_sha256'])
        require(protocol['schema']==CLASSIFICATION and protocol['parent_sha256']==PARENT_SHA and protocol['dt_s']==DT and
            protocol['numerical_comparison_limits']==LIMITS and manifest['source_sha256']==protocol['source_sha256'],
            'fixed parent/source/numerical contract differs')
        require(manifest['accepted'] is False and manifest['restart_allowed'] is False and manifest['accepted_steps']==0 and
            manifest['accepted_time_advanced_s']==0. and not list(d.rglob('accepted.npz')) and
            read(self.bind(d/'accepted_history.json'))==[],'diagnostic improperly accepted or advanced process')
        p=self.bind(protocol['parent_checkpoint'],PARENT_SHA);parent,state=checkpoint(p)
        require(state['accepted'] is True and state['restart_allowed'] is True and state['step']==14 and
            abs(state['time_s']-62.5e-6)<1e-15 and state['implementation_source_sha256']==protocol['parent_implementation_source_sha256'],
            'real parent flags/old source/time differ')
        validate_interval_energy_state(state['interval_energy']);require(len(state['interval_energy']['records'])==14,'incomplete old prefix')
        initial=self.initial(parent,state)
        rows=read(self.bind(d/'candidate_history.json'));require(len(rows)==2 and [r['route'] for r in rows]==['original','inverse'],
            'two complete original/inverse candidates required for numerical qualification')
        items={route:self.candidate(route,row,parent,state,protocol) for route,row in zip(('original','inverse'),rows)}
        declared=read(self.bind(d/'numerical_screen.json'));screen=self.compare(items['original'],items['inverse'],declared)
        require((manifest['status']==DRIVER_PASS)==screen['passed'] and receipt['exit_code']==(0 if screen['passed'] else 2),
            'terminal driver status inconsistent with independently recalculated fields')
        per=[r['wall_seconds'] for r in rows]
        require(all(np.isfinite(v) and v>0 for v in per),'invalid complete candidate costs')
        return dict(status='PASS_SAVED_INVERSE_AB_NUMERICAL_RECHECK_NOT_PROCESS_ACCEPTANCE' if screen['passed'] else
            'PASS_READBACK_OF_FAILED_AB_NOT_QUALIFIED',numerical_qualification_passed=screen['passed'],
            parent_checkpoint=rel(p),parent_sha256=PARENT_SHA,parent_time_s=state['time_s'],parent_global_step=14,
            candidate_time_s=protocol['candidate_time_s'],initial_linear_recalculation=initial,numerical_screen=screen,
            per_route={route:dict(checkpoint_sha256=items[route]['row']['sha256'],
                force_residual_N=items[route]['force_residual_N'],force_limit_N=items[route]['force_limit_N'],
                energy_prefix_records=15,iteration_log_audit=items[route]['logs']) for route in items},
            cost=dict(complete_original_candidate_s=per[0],complete_inverse_candidate_s=per[1],
                actual_original_over_inverse_ratio=per[0]/per[1],stage_wall_seconds=manifest['wall_seconds'],
                external_wall_seconds=receipt['wall_seconds'],registered_forecast_s=protocol['cost']['forecast_s'],
                comment='One ordered A/B with I/O and full gates. Stored iterate map times exclude some line-search maps; no precise decomposition or production-speed claim.'),
            limitations=['Only saved initial A/H/rhs permit independent matrix backward recomputation; later direction logs do not preserve A/H.',
                'No tangent, full force, thermal/material map or geometry certificate is rerun; those artifacts are hash/identity bound.',
                'Complete-map leaf summaries contain more fields than the saved full journal; those additional values cannot all be independently reconstructed.',
                'This same-dt algorithm check does not test temporal convergence, deposition or experiment agreement.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();directory=project(args.input);out=project(args.output)
    require(out.is_relative_to(ROOT/'outputs') and out!=ROOT/'outputs' and out!=directory and not out.is_relative_to(directory),
        'reader output must be a separate new outputs directory')
    out.mkdir(parents=True,exist_ok=False);reader=Reader(directory);start=time.perf_counter();code=2
    report=dict(status='FAIL_SAVED_INVERSE_AB_READBACK',accepted=False,restart_allowed=False,accepted_steps=0,
        accepted_time_advanced_s=0.,physical_solver_evaluations=0,nnps_evaluations=0,geometry_evaluations=0)
    try:
        reader.bind(Path(__file__).resolve());reader.bind(ROOT/'src/afsd_backend/pysph/interval_acceptance.py')
        report.update(reader.run());code=0
    except Exception as error:report.update(error=repr(error),traceback=traceback.format_exc())
    finally:
        for name,digest in reader.evidence.items():
            if sha(ROOT/name)!=digest:report.update(status='FAIL_SAVED_INVERSE_AB_READBACK',error='input changed during audit');code=2
        report.update(exit_code=code,wall_seconds=time.perf_counter()-start,evidence_sha256=reader.evidence,
            command=sys.argv,input_directory=rel(directory))
        write(out/'summary.json',report)
    print(json.dumps(dict(status=report['status'],exit_code=code,output=rel(out),wall_seconds=report['wall_seconds'])));return code


if __name__=='__main__':raise SystemExit(main())
