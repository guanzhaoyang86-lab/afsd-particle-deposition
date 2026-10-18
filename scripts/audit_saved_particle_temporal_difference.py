"""Read one particle's stored temporal-screen history; never evaluate physics."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np


ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',required=True);parser.add_argument('--particle-id',type=int,required=True)
    parser.add_argument('--pair',type=int,required=True);parser.add_argument('--output',required=True)
    options=parser.parse_args();directory=Path(options.directory).resolve();output=Path(options.output).resolve()
    if not directory.is_relative_to(ROOT/'outputs') or not output.is_relative_to(ROOT/'outputs') or output.exists():
        raise ValueError('repository output paths and unused destination required')
    bound={Path(__file__).resolve().relative_to(ROOT).as_posix():sha(Path(__file__))}
    def bind(path):
        name=path.relative_to(ROOT).as_posix();digest=sha(path)
        if name in bound and bound[name]!=digest:raise ValueError('input drift')
        bound[name]=digest;return path
    def read(path):return json.loads(bind(path).read_text(encoding='utf-8'))
    history=[]
    for number in range(1,options.pair+1):
        screen=read(directory/'fine'/f'pair_{number:03d}_temporal_screen.json')
        row=next(v for v in screen['local_screen']['rows'] if v['particle_id']==options.particle_id)
        history.append(dict(pair=number,time_s=screen['time_s'],endpoint_graphs_match=screen['endpoint_graphs_match'],**row))
    cases=[]
    choices=[('fine',2*options.pair-2),('fine',2*options.pair-1),('fine',2*options.pair),('coarse',options.pair)]
    for stage,number in choices:
        folder=directory/stage/f'step_{number:03d}';path=bind(folder/'NOT_RESTARTABLE_candidate.npz')
        with np.load(path,allow_pickle=True) as data:
            state=data['solver_data'].item();arrays=data['particles'].item()['solid']['arrays']
        if state['restart_allowed'] is not False:raise ValueError('expected explicit raw candidate')
        index=int(np.flatnonzero(arrays['particle_id']==options.particle_id)[0]);stress=arrays['S'].reshape(-1,3,3)[index]
        contact=read(folder/'FULL_UNCOMMITTED_candidate.json')['contact']
        with np.load(bind(folder/'visible_graph.npz'),allow_pickle=False) as graph:
            ids=graph['particle_id'];i,j=graph['i'],graph['j'];mask=(ids[i]==options.particle_id)|(ids[j]==options.particle_id)
            neighbors=np.unique(np.r_[ids[i[mask]],ids[j[mask]]]).tolist();blocked=len(graph['blocked_raw_edge_indices'])
        cases.append(dict(stage=stage,step=number,time_s=state['time_s'],accepted_checkpoint_present=(folder/'accepted.npz').exists(),
            raw_candidate_checkpoint=path.relative_to(ROOT).as_posix(),actuator=int(arrays['actuator'][index]),free=state['free'][index].tolist(),
            particle_fields={key:float(arrays[key][index]) for key in ('x','y','z','temperature','ep','ep_rate','strain_rate','plastic_work','return_dissipation')},
            stress_Pa=stress.tolist(),stress_frobenius_norm_Pa=float(np.linalg.norm(stress)),
            contact_active=contact['active'][index],contact_sliding=contact['sliding'][index],contact_features=contact['features'],
            contact_force_N=contact['force_N'][index],contact_signed_energy_J={key:contact[key][index] for key in
                ('friction_dissipation_J','algorithmic_dissipation_J','release_energy_J')},
            active_neighbor_ids_including_self=neighbors,global_blocked_edge_count=blocked))
    first=np.asarray(cases[-2]['stress_Pa']);last=np.asarray(cases[-1]['stress_Pa'])
    cosine=float(np.sum(first*last)/(np.linalg.norm(first)*np.linalg.norm(last)))
    if not -1<=cosine<=1:raise ValueError('invalid tensor cosine; no clipping')
    result=dict(status='PASS_SAVED_PARTICLE_TEMPORAL_DIAGNOSTIC_NOT_CAUSAL_OR_PHYSICAL_VALIDATION',
        command=[sys.executable,'-B',str(Path(__file__).resolve()),*sys.argv[1:]],physical_solver_evaluations=0,nnps_evaluations=0,
        source_and_input_sha256=bound,particle_id=options.particle_id,comparison_history=history,cases=cases,
        final_stress_difference_norm_Pa=float(np.linalg.norm(first-last)),
        final_stress_norm_relative_difference=float(abs(np.linalg.norm(first)-np.linalg.norm(last))/np.linalg.norm(first)),
        final_stress_tensor_component_vector_angle_degrees=float(np.arccos(cosine)*180/np.pi),
        interpretation='Tensor component-vector angle is a diagnostic, not a material rigid rotation or a causal attribution.',
        original_failed_temporal_screen_preserved=True,accepted_time_advanced_s=0.,exit_code=0)
    if any(sha(ROOT/name)!=value for name,value in bound.items()):raise ValueError('input changed during audit')
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','final_stress_difference_norm_Pa','final_stress_norm_relative_difference','final_stress_tensor_component_vector_angle_degrees')},indent=2))


if __name__=='__main__':main()
