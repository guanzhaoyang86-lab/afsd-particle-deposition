"""Read-only independent inventory/status checks; not a solver reference."""
from pathlib import Path
import gzip
import hashlib
import json
import re

ROOT=Path(__file__).resolve().parents[1]
LOGS=ROOT/'logs/process_event_20261003'
OUT=ROOT/'outputs/process_event_20261003'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def require(value,message):
    if not value:raise ValueError(message)


def validate(directory,report):
    for rel,h in dict(report.get('input_sha256',{}),**report['source_sha256']).items():
        p=(ROOT/rel).resolve();require(p.is_relative_to(ROOT) and sha(p)==h,'source/input drift '+rel)
    for rel,h in report['artifact_sha256'].items():
        p=(directory/rel).resolve();require(p.is_relative_to(directory) and sha(p)==h,'artifact drift '+rel)


def main():
    pre=read(LOGS/'preflight_refined/report.json')
    stable=read(LOGS/'preflight_stable/report.json')
    run=read(OUT/'response_probe/manifest.json')
    for directory,report in ((LOGS/'preflight_refined',pre),(LOGS/'preflight_stable',stable),(OUT/'response_probe',run)):
        validate(directory,report)
    require(pre['status']=='FAIL_PREFLIGHT' and pre['tests_run']==146 and
        pre['tests_passing']==144 and pre['tests_failed']==2 and pre['tests_errors']==0,'unexpected preflight')
    failed=[r for r in pre['tests'] if r['tests_failed']]
    require([r['name'] for r in failed]==['event_response','event_response_refined'],'different failed tests')
    require(stable['tests_run']==139 and stable['tests_failed']==1,'old gate drift')
    wiring=(LOGS/'event_refined_wiring_initial.log').read_text(encoding='utf-8')
    require(re.search(r'^Ran 1 test in ',wiring,re.M) and re.search(r'^OK\s*$',wiring,re.M),'isolated wiring not PASS')
    require(run['status']=='BLOCKED_NO_NEW_ACCEPTED_STATE' and run['accepted_steps']==0 and
        run['total_evaluations']==0 and run['accepted_time_s']==run['source_time_s'],'false advance')
    require(read(OUT/'response_probe/accepted_history.json')==[],'unaccepted probe published history')
    parent=run['parent_reference'];require(sha(ROOT/parent['path'])==parent['sha256'],'parent drift')
    with gzip.open(OUT/'response_probe/particle_000001.json.gz','rt',encoding='utf-8') as f:r=json.load(f)
    require(not r['available'] and not r['accepted'] and r['B'] is None and r['C'] is None,'invalid derivative availability')
    require(r['particle_id']==1 and r['kick_calls']==5 and r['reason']=='ONE_SIDED_DIFFERENCES_INCONSISTENT','probe failure mismatch')
    c=r['column_evidence'][0]
    # Independent subtraction of saved derivatives; never call a tested response.
    gap=abs(c['forward_half_h_s_kg'][1][1]-c['backward_half_h_s_kg'][1][1])
    bound=c['one_sided_allowed_s_kg'][1][1]
    require(gap>bound and c['refinement_pass'] and not c['one_sided_pass'],'failure not supported')
    for name in ('response_refined','cost_probe','root_attempt'):
        require(not (OUT/name).exists(),'unexpected later run '+name)
    for path in OUT.rglob('accepted.npz'):raise ValueError('unexpected physical state '+str(path))
    print(json.dumps(dict(status='PASS_EVIDENCE_READBACK_PROJECT_INCOMPLETE',
        final_numerical_gate='FAIL_PREFLIGHT',distinct_tests_executed=147,passing=145,failed=2,
        isolated_wiring='PASS_NOT_GATE_OVERRIDE',new_physical_evaluations=0,new_accepted_steps=0,
        accepted_time_s=run['accepted_time_s'],parent_reference=parent,
        first_unavailable=dict(particle_id=1,contact_id='1:0',gap_Cyx_s_kg=gap,bound_s_kg=bound),
        refined_cloud='NOT_RUN_FAILED_REGRESSION',cost_probe='NOT_RUN',fullroot='NOT_RUN',
        time_sensitivity='NOT_RUN',space_sensitivity='NOT_RUN',production_authorized=False,
        source_artifact_and_parent_hashes='PASS',old_failure_retained=True,project_complete=False),indent=2))


if __name__=='__main__':main()
