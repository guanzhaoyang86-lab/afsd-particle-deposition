"""Fresh original acceptance gates for an independent expanded candidate only.

This helper neither publishes nor marks its detached result restartable. The
process controller must still complete geometry and same-parent pair screens.
"""
from pathlib import Path
import numpy as np
import expanded_hot_domain as domain
from afsd_backend.pysph.discrete_process import accepted_copy as original_accepted_copy
from afsd_backend.pysph.interval_acceptance import validate_interval_energy_state, validate_bore_crossing_state

QUALIFICATION = 'EXPANDED_8MM_INITIAL_ELASTIC_RANK_AND_SAME_MAP_V1'
PASS = 'PASS_EXPANDED_INITIAL_ELASTIC_RANK_AND_SAME_MAP_NOT_PHYSICAL_VALIDATION'


def accepted_copy(trial, result, qualification):
    """Return (pa,state,fresh) after original fresh-force/raw/all-prefix gates.

    qualification is a path to an actual completed independent receipt, not an
    assertion dictionary. Neither old spatial transitions nor old parents are
    fabricated. Coarse states may be chained diagnostically but stay unpublished.
    """
    from qualify_expanded_hot_domain import read, sha, verify, immutable_hashes, INITIAL_SHA, DT
    require=domain.require;path=Path(qualification)
    if not path.is_absolute():path=domain.ROOT/path
    path=path.resolve();require(path.is_relative_to(domain.ROOT),'qualification outside project')
    q=read(path);supervision=read(path.parent/'supervision.json')
    require(q.get('status')==PASS and q.get('qualification')==QUALIFICATION and q.get('qualification_passed') is True,
        'independent initial rank and map qualification missing')
    require(q.get('accepted') is False and q.get('restart_allowed') is False and q['accepted_time_advanced_s']==0.,
        'initial qualification must not claim process acceptance')
    require(supervision['exit_code']==0 and supervision['before_after_hashes_match'] and
        supervision['summary_sha256']==sha(path) and supervision['log_sha256']==sha(path.parent/'run.log'),
        'qualification external receipt/hash mismatch')
    require(q['initialization_sha256']==INITIAL_SHA and q['reference_dt_s']==DT and q['particle_count']==918 and
        q['free_displacement_dof']==1961,'different independent model qualification')
    require(q['initial_elastic_rank']['status']=='PASS_INITIAL_FREE_ELASTIC_RANK_NOT_FUTURE_FLOW_STABILITY' and
        q['initial_same_map']['status']=='PASS_INDEPENDENT_INITIAL_SAME_GRAPH_COMPLETE_MAP_AND_TANGENT',
        'incomplete independent rank or same-map receipt')
    verify(q['source_and_input_sha256']);verify(q['output_sha256'])
    state=trial.state
    require(state.get('classification')==domain.CLASSIFICATION and state.get('model_scope')==domain.SCOPE and
        state.get('operator_variant')=='visible' and state.get('independent_initialization') is True,
        'independent expanded-visible model required')
    require('spatial_model_transition' not in state and state.get('branch_origin_sha256',INITIAL_SHA)==INITIAL_SHA,
        'independent candidate cannot use a fabricated old transition or parent')
    if 'continuation_parent_sha256' in state:
        parent=state['continuation_parent_sha256']
        require(isinstance(parent,str) and len(parent)==64 and all(c in '0123456789abcdef' for c in parent),
                'invalid direct continuation parent hash')
        require(state.get('window_origin_sha256')==INITIAL_SHA,
                'startup continuation must keep its real independent window origin')
    require(state['parameters']==q['parameters'] and state['parameter_provenance']==q['parameter_provenance'] and
        domain.clean(state['independent_domain_initialization'])==q['independent_domain_initialization'],
        'independent model provenance differs from qualification')
    require(immutable_hashes(trial.pa,state)==q['immutable_array_sha256'],'independent identity/axes/mass changed')
    validate_interval_energy_state(state['interval_energy']);validate_bore_crossing_state(state['bore_observation'])
    require(state['interval_energy']['steps']==state['step'],'complete accepted/diagnostic energy prefix missing')
    if state['step']==0:
        require(state['time_s']==0. and state.get('accepted') is False and state.get('restart_allowed') is False and
                'branch_origin_sha256' not in state and 'continuation_parent_sha256' not in state,
                'independent time-zero state cannot already be accepted')
        domain.initial_geometry_and_field_audit(trial.pa,state)
    else:require(state['time_s']>0. and state.get('expanded_initial_qualification',{}).get('receipt_sha256')==sha(path),
        'continued independent candidate lacks its initial qualification binding')
    nextpa,nextstate,fresh=original_accepted_copy(trial,result)
    require(immutable_hashes(nextpa,nextstate)==q['immutable_array_sha256'],'fresh acceptance changed immutable fields')
    require(nextstate['interval_energy']['records'][:-1]==state['interval_energy']['records'],
            'fresh acceptance reset the complete energy prefix')
    require(nextstate['step']==state['step']+1 and np.isclose(nextstate['time_s'],state['time_s']+result['dt_s'],rtol=0.,atol=0.),
            'fresh acceptance advanced an unexpected step/time')
    nextstate.update(accepted=False,restart_allowed=False,expanded_initial_qualification=dict(
        qualification=QUALIFICATION,receipt=path.relative_to(domain.ROOT).as_posix(),receipt_sha256=sha(path),
        initialization_sha256=INITIAL_SHA,publication_requires_same_parent_temporal_screen=True))
    verify(q['source_and_input_sha256'])
    return nextpa,nextstate,fresh
