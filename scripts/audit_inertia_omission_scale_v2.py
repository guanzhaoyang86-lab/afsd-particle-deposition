"""Qualify available force scales before interpreting inertial scale ratios.

Preserves the initial arithmetic audit verbatim. In a free group without
external contact, saved resultant internal force can itself be merely the QS
equilibrium residual. Such a group has no resolved force denominator here.
Use the saved solver's componentwise force tolerance times sqrt(group DOFs)
only to mark an unavailable scale, never as a substituted denominator. No
physical acceptance threshold, model or saved trajectory is changed.
"""
from pathlib import Path
import copy
import json
import numpy as np
import audit_inertia_omission_scale as original

ROOT = original.ROOT
OLD = original.DEFAULT_OUTPUT.with_suffix('.json')
OUT = ROOT/'outputs/coupled_cost_20261004/inertia_omission_30ms_v2.json'


def main():
    original.require(not OUT.exists(), 'refuse qualified diagnostic overwrite')
    initial = json.loads(OLD.read_text(encoding='utf-8'))
    for name,digest in dict(initial['source_sha256'],**initial['input_sha256']).items():
        original.require(original.reader.sha(ROOT/name)==digest, 'original diagnostic input drift: '+name)
    arrays = OLD.with_suffix('.npz')
    original.require(original.reader.sha(arrays)==initial['arrays_sha256'], 'raw diagnostic arrays drift')
    summary = original.reader.read(original.READBACK/'summary.json')
    with np.load(arrays,allow_pickle=False) as saved:
        free=saved['physical_free_dofs'].copy()
        fi,fc=saved['saved_internal_force_N'].copy(),saved['saved_contact_force_N'].copy()
        inertia=saved['inertial_force_proxy_N'].copy()
    with np.load(original.READBACK/'accepted_trajectory.npz',allow_pickle=False) as saved:
        component=saved['particle_property__component_id'][0,:,0].copy()
    masks=dict(all_dofs=np.ones_like(free,bool),physical_free_dofs=free,
        bar_free_dofs=free&(component==1)[:,None],substrate_free_dofs=free&(component==0)[:,None])
    report=copy.deepcopy(initial)
    for n,row in enumerate(report['rows']):
        limit=summary['rows'][n]['force_limit_N']
        for name,group in row['groups'].items():
            envelope=float(np.sqrt(masks[name].sum())*limit)
            resolved=group['constituent_force_scale_l2_N']>envelope
            group.update(registered_equilibrium_tolerance_envelope_l2_N=envelope,
                constituent_scale_resolved_above_tolerance=bool(resolved))
            if not resolved:
                group['inertial_to_constituent_ratio']=None
                group['unavailable_reason']='Saved nodal constituent scale is at/below equilibrium tolerance; do not divide by a near-zero residual.'
        particle_scale=np.sqrt(np.sum((fi[n]*free)**2+(fc[n]*free)**2,axis=1))
        particle_inertia=np.linalg.norm(inertia[n]*free,axis=1)
        cutoff=.01*float(particle_scale.max())
        local_envelope=np.sqrt(np.count_nonzero(free,axis=1))*limit
        supported=(particle_scale>local_envelope)&(particle_scale>=cutoff)&(particle_scale>0.)
        row.update(supported_particle_count=int(supported.sum()),
            unsupported_or_small_scale_particle_count=int((~supported).sum()),
            supported_pointwise_ratio_statistics=original.statistics(particle_inertia[supported]/particle_scale[supported]))
    rows=report['rows']
    report['aggregate_ratios']={name:dict(original.statistics([row['groups'][name]['inertial_to_constituent_ratio'] for row in rows]),
        unresolved_node_count=sum(row['groups'][name]['inertial_to_constituent_ratio'] is None for row in rows)) for name in masks}
    report['stage_ratios']={}
    for stage in sorted({row['stage'] for row in rows}):
        subset=[row for row in rows if row['stage']==stage]
        report['stage_ratios'][stage]=dict(nodes=len(subset),start_s=subset[0]['time_s'],end_s=subset[-1]['time_s'],
            free_ratio_statistics=original.statistics([row['groups']['physical_free_dofs']['inertial_to_constituent_ratio'] for row in subset]))
    available=[row for row in rows if row['groups']['physical_free_dofs']['inertial_to_constituent_ratio'] is not None]
    report['worst_free_ratio_node']=max(available,key=lambda row:row['groups']['physical_free_dofs']['inertial_to_constituent_ratio'])
    report['dt_transition_nodes']=[row for row in rows if row['dt_changed']]
    report['constant_dt_free_ratio_statistics']=original.statistics([row['groups']['physical_free_dofs']['inertial_to_constituent_ratio'] for row in rows if not row['dt_changed']])
    report['post_0_5ms_free_ratio_statistics']=original.statistics([row['groups']['physical_free_dofs']['inertial_to_constituent_ratio'] for row in rows if row['time_s']>=.0005])
    report['status']='COMPUTED_WITH_UNRESOLVED_SUBGROUPS_READ_ONLY_INERTIAL_SCALE_NOT_VALIDATION'
    report['force_scale_resolvability']='A scale must exceed sqrt(group DOF count)*saved componentwise force tolerance to admit a ratio; no tolerance is substituted for the physical scale.'
    report['pointwise_diagnostic_subset']='Constituent nodal scale >=1% of node maximum AND above its free-component tolerance envelope. Excluded points remain unresolved.'
    report['interpretation']='Global/bar force scale is available. Substrate free resultant internal force is not a resolved denominator: no substrate inertia-omission conclusion follows.'
    report['limitation_added']='The older v1 raw substrate ratios used a denominator near equilibrium residual. They are not interpretable and are superseded by null/unavailable values in this report.'
    report['superseded_initial_arithmetic_report']=dict(path=OLD.relative_to(ROOT).as_posix(),sha256=original.reader.sha(OLD),old_file_modified=False)
    report['raw_arrays_file']=arrays.relative_to(ROOT).as_posix()
    report['source_sha256'][Path(__file__).relative_to(ROOT).as_posix()]=original.reader.sha(Path(__file__))
    report['input_sha256'][OLD.relative_to(ROOT).as_posix()]=original.reader.sha(OLD)
    report['input_sha256'][arrays.relative_to(ROOT).as_posix()]=original.reader.sha(arrays)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(output=str(OUT),aggregate=report['aggregate_ratios'],post_0_5ms=report['post_0_5ms_free_ratio_statistics'])))


if __name__=='__main__':
    main()
