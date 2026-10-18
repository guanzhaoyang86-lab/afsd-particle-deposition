"""One bounded saved342 material-P action; zero physical evaluations."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,re,sys,threading,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.material_stage_response import MaterialStageResponse
FIELDS=('position_m','velocity_m_s','density_kg_m3','temperature_K','mean_velocity_m_s')
DATA=ROOT/'outputs/process_bridge_20261003/real342_initial'
OLD=ROOT/'logs/process_bridge_action_repair1_20261004'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def gz(p):
    with gzip.open(p,'rt',encoding='utf-8') as f:return json.load(f)
def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {str(k):plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [plain(x) for x in v]
    return v
def write(p,v):Path(p).write_text(json.dumps(plain(v),indent=2,allow_nan=False)+'\n',encoding='utf-8')
def require(v,s):
    if not v:raise ValueError(s)
def bind(values,base,watched):
    for name,digest in values.items():
        path=(base/name).resolve()
        require(path.is_relative_to(base.resolve()) and sha(path)==digest,'evidence/source drift '+name)
        watched[path.relative_to(ROOT).as_posix()]=digest

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--material-tests',type=Path,required=True)
    ap.add_argument('--schur-tests',type=Path,required=True)
    args=ap.parse_args();out=args.output.resolve()
    require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);started=time.perf_counter()
    def timeout():
        write(out/'TIMEOUT.json',dict(status='INCOMPLETE55S',accepted=False,new_evaluations=0));os._exit(124)
    timer=threading.Timer(55,timeout);timer.daemon=True;timer.start()
    report=dict(status='RUNNING',accepted=False,restart_allowed=False,can_publish=False,
        new_evaluations=0,new_afsd_steps=0,new_transport_events=0,production_authorized=False)
    watched={};code=2
    try:
        previous=read(ROOT/'docs/PROCESS_BRIDGE_FREEZE.json')
        require(len(previous['files_sha256'])==4745,'unexpected parent freeze')
        for name,digest in previous['files_sha256'].items():
            require(sha(ROOT/name)==digest or name in ALLOWED,'prior evidence changed '+name)
        watched['docs/PROCESS_BRIDGE_FREEZE.json']=sha(ROOT/'docs/PROCESS_BRIDGE_FREEZE.json')
        test_records=[]
        for arg,required in ((args.material_tests,'src/afsd_backend/pysph/material_response.py'),
                             (args.schur_tests,'src/afsd_backend/pysph/material_stage_response.py')):
            path=arg.resolve();require(path.is_relative_to(ROOT),'local receipt required');r=read(path)
            require(r['return_code']==0 and r['sources_unchanged'] and
                    r['source_hashes_before']==r['source_hashes_after'],'component tests failed/source drift')
            require(required in r['source_hashes_before'],'component not bound by tests')
            bind(r['source_hashes_before'],ROOT,watched)
            logs=[p for p in path.parent.glob('*.log') if sha(p)==r['log_sha256']]
            require(len(logs)==1,'unique bound component log required')
            text=logs[0].read_text(encoding='utf-8');counts=re.findall(r'^Ran (\d+) tests? in ',text,re.M)
            require(len(counts)==1 and re.search(r'^OK\s*$',text,re.M) is not None,'test outcome incomplete')
            watched[path.relative_to(ROOT).as_posix()]=sha(path)
            watched[logs[0].relative_to(ROOT).as_posix()]=sha(logs[0])
            test_records.append(dict(path=path.relative_to(ROOT).as_posix(),tests=int(counts[0])))
        actual=read(DATA/'report.json');old=read(OLD/'report.json')
        require(actual['status']=='NONSMOOTH_NOT_CERTIFIED' and actual['actual_evaluations']==4,'wrong real342 archive')
        for folder,r in ((DATA,actual),(OLD,old)):
            watched[(folder/'report.json').relative_to(ROOT).as_posix()]=sha(folder/'report.json')
            bind(r['source_sha256'],ROOT,watched);bind(r['artifact_sha256'],folder,watched)
        for name in ('scripts/probe_material_response_action.py','scripts/run_material_response_action.py',
                     'docs/MATERIAL_RESPONSE_ACTION_CONTRACT.md',
                     'src/afsd_backend/pysph/material_stage_response.py'):
            watched[name]=sha(ROOT/name)
        saved={label:gz(DATA/(label+'_candidate.json.gz')) for label in ('base','plus','minus')}
        parent=gz(DATA/'parent_state.json.gz');p=read(DATA/'configuration.json')['parameters']
        before=sha(DATA/'parent_state.json.gz')
        ids=np.asarray(parent['properties']['particle_id']);n=len(ids)
        tick=time.perf_counter()
        pre=MaterialStageResponse.from_journal(parent['properties'],parent['state'],p,
                                               saved['base']['candidate'],saved['base']['evaluation'])
        build=time.perf_counter()-tick;applied={};cost={}
        for label in ('plus','minus','base'):
            tick=time.perf_counter()
            applied[label]=pre.apply(saved[label]['candidate'],saved[label]['evaluation']['map'])
            cost[label]=time.perf_counter()-tick
        npz=np.load(OLD/'action_arrays.npz',allow_pickle=False)
        require(np.array_equal(ids,npz['particle_id']),'old comparison IDs differ')
        dc={k:np.asarray(saved['plus']['candidate'][k])-saved['minus']['candidate'][k] for k in FIELDS}
        traces=saved['base']['evaluation']['boundary']['traces']
        nc=np.asarray([len(t['contact_rows']) for t in traces]);zero=np.asarray([bool(t['result']['events']) for t in traces])
        ep=np.asarray(saved['base']['evaluation']['material']['ep'])-parent['properties']['ep']
        groups=dict(all=np.ones(n,bool),plastic=ep>0,no_contact=nc==0,multiple_contact=nc>1,zero_event=zero)
        scales=dict(position_m=1e-12,velocity_m_s=1e-9,mean_velocity_m_s=1e-9,
                    density_kg_m3=np.asarray(parent['properties']['rho'])*1e-10,temperature_K=5e-12)
        arrays={'particle_id':ids};comparison={}
        for key in FIELDS:
            corr=applied['plus']['correction'][key]-applied['minus']['correction'][key]
            error=dc[key]+corr
            arrays['delta_candidate_'+key]=dc[key];arrays['delta_correction_'+key]=corr;arrays['action_error_'+key]=error
            arrays['delta_residual_'+key]=np.asarray(saved['plus']['evaluation']['residual'][key])-saved['minus']['evaluation']['residual'][key]
            require(np.array_equal(dc[key],npz['delta_candidate_'+key]),'changed finite direction')
            comparison[key]={}
            for group,mask in groups.items():
                row={}
                for name,value in (('input',dc[key]),('new_material',error),
                                   ('old_structural',npz['structural_action_error_'+key]),
                                   ('old_continuity',npz['continuity_action_error_'+key])):
                    v=value[mask];scale=np.broadcast_to(scales[key],value.shape)[mask]
                    row[name]=dict(inf=float(np.max(abs(v),initial=0)),l2=float(np.linalg.norm(v)),
                        scaled_inf=float(np.max(abs(v/scale),initial=0)))
                den=row['input']['l2']
                row['new_to_input_l2']=row['new_material']['l2']/den if den else None
                row['new_to_old_continuity_l2']=(row['new_material']['l2']/row['old_continuity']['l2']
                                                if row['old_continuity']['l2'] else None)
                comparison[key][group]=row
        np.savez_compressed(out/'action_arrays.npz',**arrays)
        # An algebraic candidate is never an evaluated or publishable endpoint.
        with gzip.open(out/'UNEVALUATED_BASE_PROPOSAL.json.gz','wt',encoding='utf-8') as f:
            json.dump(plain(dict(accepted=False,restart_allowed=False,can_publish=False,
                classification='UNEVALUATED_NUMERICAL_PROPOSAL_NOT_PHYSICAL_STATE',
                result=applied['base'])),f,allow_nan=False)
        for name,digest in watched.items():require(sha(ROOT/name)==digest,'watched file changed '+name)
        require(sha(DATA/'parent_state.json.gz')==before,'parent changed')
        report.update(status='COMPLETE_SAVED_ACTION_NOT_ROOT_OR_PROCESS',source_sha256=watched,
            component_tests=test_records,particle_count=n,parent_time_s=parent['state']['t'],dt_s=actual['physical_dt_s'],
            new_build_seconds=build,apply_seconds=cost,comparison=comparison,
            backward_errors={k:v['diagnostics']['backward_error_inf'] for k,v in applied.items()},
            matrix_asymmetry_inf_kg=applied['base']['diagnostics']['matrix_asymmetry_inf_kg'],
            model_metadata=pre.metadata,parent_unchanged=True,
            source_classification=actual['status'],physical_fields_updated=False)
        code=0
    except Exception as ex:report.update(status='FAIL_SAVED_ACTION',error=repr(ex))
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started,
        artifact_sha256={p.name:sha(p) for p in out.glob('*') if p.is_file() and p.name!='report.json'})
    write(out/'report.json',report)
    print(json.dumps(plain({k:v for k,v in report.items() if k not in ('source_sha256','comparison','model_metadata')})),flush=True)
    return code

if __name__=='__main__':raise SystemExit(main())

