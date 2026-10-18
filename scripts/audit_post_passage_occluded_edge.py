"""Independent Decimal readback of the rejected 39.875 ms NNPS edge.

Only stored checkpoints, journals, edges and a prescribed pose are read.
No physical residual, NNPS, gradient, force, contact or thermal update occurs.
The strict segment/solid intersections are recomputed with 70-digit Decimal
arithmetic on the saved world coordinates and double-precision prescribed pose.
"""
from pathlib import Path
from decimal import Decimal, localcontext
import json
import numpy as np
import readback_discrete_process as reader

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/coupled_cost_20261004/batched/post_passage_short_response/coarse'
OUT=ROOT/'outputs/coupled_cost_20261004/post_passage_occluded_edge_audit.json'
D=Decimal.from_float


def local_decimal(position,origin,rotation):
    return [[sum((D(float(point[j]))-D(float(origin[j])))*D(float(rotation[j,k])) for j in range(3))
             for k in range(3)] for point in position]


def segment_intersection(a,b,radius,half_bore,height,inset=0.):
    zero,one=Decimal(0),Decimal(1)
    delta=D(float(inset));r=D(float(radius))-delta;q=D(float(half_bore))+delta
    bottom,top=delta,D(float(height))-delta
    direction=[right-left for left,right in zip(a,b)]
    roots=[(zero,'segment_start'),(one,'segment_end')]
    plane_roots={}
    for axis,value,name in ((0,-q,'BORE_X_MINUS'),(0,q,'BORE_X_PLUS'),
            (1,-q,'BORE_Y_MINUS'),(1,q,'BORE_Y_PLUS'),(2,bottom,'SHOULDER'),(2,top,'TOP')):
        if direction[axis]:
            root=(value-a[axis])/direction[axis];plane_roots[name]=root
            if zero<root<one:roots.append((root,name))
    aa=sum(v*v for v in direction[:2]);bb=2*sum(a[k]*direction[k] for k in range(2));cc=sum(v*v for v in a[:2])-r*r
    disc=bb*bb-4*aa*cc
    if aa>0 and disc>=0:
        for sign in (-1,1):
            root=(-bb+sign*disc.sqrt())/(2*aa)
            if zero<root<one:roots.append((root,'OUTER_CYLINDER'))
    cuts=sorted(set(value for value,_ in roots));inside=[]
    for lo,hi in zip(cuts,cuts[1:]):
        mid=(lo+hi)/2;point=[a[k]+mid*direction[k] for k in range(3)]
        if bottom<point[2]<top and sum(v*v for v in point[:2])<r*r and max(abs(v) for v in point[:2])>q:
            if inside and inside[-1][1]==lo:inside[-1][1]=hi
            else:inside.append([lo,hi])
    length=sum(v*v for v in direction).sqrt()
    rows=[]
    for lo,hi in inside:
        mid=(lo+hi)/2;point=[a[k]+mid*direction[k] for k in range(3)]
        rows.append(dict(start_fraction=str(lo),end_fraction=str(hi),length_m=str((hi-lo)*length),
            start_features=sorted(name for root,name in roots if root==lo),
            end_features=sorted(name for root,name in roots if root==hi),
            midpoint_local_m=[str(v) for v in point],midpoint_above_bottom_m=str(point[2]-bottom),
            midpoint_outside_square_m=str(max(abs(v) for v in point[:2])-q),
            midpoint_inside_cylinder_radial_m=str(r-sum(v*v for v in point[:2]).sqrt())))
    zroot=plane_roots['SHOULDER'];yroot=plane_roots['BORE_Y_PLUS']
    y_at_z=a[1]+zroot*direction[1]-q
    return dict(precision_decimal_digits=70,geometric_inset_m=inset,edge_length_m=str(length),
        strict_inside_intervals=rows,maximum_strict_length_m=str(max([(hi-lo)*length for lo,hi in inside],default=zero)),
        shoulder_fraction=str(zroot),bore_y_plus_fraction=str(yroot),
        y_minus_positive_bore_wall_at_shoulder_intersection_m=str(y_at_z),
        shoulder_to_bore_root_separation=str(yroot-zroot))


def main():
    reader.require(not OUT.exists(),'refuse failure audit overwrite')
    inputs={}
    def bind(path):
        inputs[path.relative_to(ROOT).as_posix()]=reader.sha(path);return path
    def read(path):return reader.read(bind(path))
    def load(path):
        data=reader.load(str(bind(path)));return data['arrays']['solid'],data['solver_data']
    meta=read(BASE/'manifest.json');history=read(BASE/'candidate_history.json')
    reader.require(meta['status']=='FAILED_NO_NEW_ACCEPTED_STATE' and meta['accepted_steps']==0 and
        len(history)==meta['candidate_steps']==19 and not list(BASE.glob('step_*/accepted.npz')),
        'rejected post-passage branch must remain unaccepted')
    reader.require(not (BASE.parent/'fine').exists(),'fine must not exist after failed coarse')
    failure=read(BASE/'step_020/solid_edge_visibility.json')
    reader.require(failure['status']=='FAIL_SOLID_OCCLUDED_ENDPOINT_EDGE' and failure['occluded_edge_count']==1,
        'expected single-edge actual rejection absent')
    failed_edge=failure['occluded_edges'][0]
    pids=[failed_edge['particle_id_i'],failed_edge['particle_id_j']]
    snapshots=[('accepted_30ms',ROOT/'outputs/coupled_cost_20261004/batched/finite_shoulder_refinement/fine/step_128/accepted.npz'),
        ('accepted_first_passage',ROOT/meta['parent_checkpoint']),
        ('coarse_previous_39_75ms',BASE/'step_019/NOT_RESTARTABLE_candidate.npz'),
        ('coarse_rejected_39_875ms',BASE/'step_020/NOT_RESTARTABLE_candidate.npz')]
    results=[]
    with localcontext() as context:
        context.prec=70
        for label,path in snapshots:
            pa,state=load(path);folder=path.parent;ids=pa.particle_id
            index=[int(np.flatnonzero(ids==pid)[0]) for pid in pids]
            positions=np.column_stack((pa.x,pa.y,pa.z));p=state['parameters'];pose=reader.expected_pose(p,state['time_s'])
            local=local_decimal(positions[index],pose['origin'],pose['rotation'])
            evidence_path=bind(folder/'candidate_evidence.npz')
            with np.load(evidence_path,allow_pickle=False) as evidence:
                reader.require(np.array_equal(ids,evidence['particle_id']),'saved NNPS edge particle order mismatch')
                matches=np.flatnonzero((evidence['i']==index[0])&(evidence['j']==index[1]))
                edges=len(evidence['i'])
            reader.require(len(matches)==1,'audited edge absent or duplicated in an actual saved graph')
            strict=segment_intersection(local[0],local[1],p['tool_radius'],p['bore_width']/2,p['tool_top']-p['layer'])
            inset=segment_intersection(local[0],local[1],p['tool_radius'],p['bore_width']/2,p['tool_top']-p['layer'],1e-10)
            visibility=read(folder/'solid_edge_visibility.json')
            contact=read(folder/'FULL_UNCOMMITTED_candidate.json')['contact']
            particles=[]
            for pid,k,point in zip(pids,index,local):
                particles.append(dict(particle_id=pid,raw_index=k,component_id=int(pa.component_id[k]),
                    material_id=int(pa.material_id[k]),tracer_type=int(pa.tracer_type[k]),wholebar_tag=int(pa.wholebar_tag[k]),
                    world_position_m=positions[k].tolist(),tool_local_position_m=[str(v) for v in point],
                    active_contact_features=[name for name,active in zip(contact['features'],contact['active'][k]) if active],
                    contact_gap_m=dict(zip(contact['features'],contact['gaps_m'][k])),
                    particle_radius_m=p['dx']/2,smoothing_length_m=float(pa.h[k]),mass_kg=float(pa.m[k])))
            item=dict(label=label,time_s=state['time_s'],global_step=state['step'],
                accepted=path.name=='accepted.npz',restart_allowed=state.get('restart_allowed',True),
                checkpoint_sha256=reader.sha(path),particles=particles,edge_count=edges,
                edge_index=int(matches[0]),edge_present_in_saved_nnps_graph=True,
                strict_intersection=strict,inset_intersection=inset,
                saved_visibility_status=visibility['status'],saved_occluded_count=visibility['occluded_edge_count'])
            results.append(item)
        final=results[-1]
        computed=float(final['strict_intersection']['maximum_strict_length_m'])
        inset_computed=float(final['inset_intersection']['maximum_strict_length_m'])
        reader.require(abs(computed-failed_edge['maximum_contiguous_inside_length_m'])<1e-16 and
            abs(inset_computed-failed_edge['inset_maximum_contiguous_inside_length_m'])<1e-16,
            'independent high-precision intersections disagree with saved guard')
        reader.require(computed>100*failure['minimum_occluded_length_m'] and inset_computed>100*failure['minimum_occluded_length_m'],
            'claimed robust occlusion is not above threshold by two orders')
        reader.require(all(not item['strict_intersection']['strict_inside_intervals'] for item in results[:-1]),
            'previous sampled states already contain this occlusion')
        midpoint=final['strict_intersection']['strict_inside_intervals'][0]
        feature_set=midpoint['start_features']+midpoint['end_features']
        reader.require(feature_set==['SHOULDER','BORE_Y_PLUS'],'different finite tool feature intersection')
    geometry=read(BASE/'step_020/geometry_audit.json');shell=read(BASE/'step_020/continuous_shell_audit.json')
    scalar=read(BASE/'step_020/scalar_batch_equivalence.json')
    report=dict(status='CONFIRMED_REAL_TOOL_SOLID_OCCLUDED_NNPS_EDGE_REJECTION_PRESERVED',
        physical_solver_evaluations=0,nnps_evaluations=0,material_updates=0,force_evaluations=0,
        independent_geometry_recomputation='70-digit Decimal segment intersections on saved positions and prescribed double pose',
        accepted_time_unchanged_s=meta['accepted_time_s'],coarse_failed_candidate_time_s=failure['endpoint_time_s'],
        accepted_checkpoint_sha256=meta['parent_sha256'],pids=pids,states=results,
        guard_comparison=dict(saved_strict_length_m=failed_edge['maximum_contiguous_inside_length_m'],
            independent_strict_length_m=computed,saved_inset_length_m=failed_edge['inset_maximum_contiguous_inside_length_m'],
            independent_inset_length_m=inset_computed,strict_length_over_guard_threshold=computed/failure['minimum_occluded_length_m'],
            inset_length_over_guard_threshold=inset_computed/failure['minimum_occluded_length_m'],
            strict_guard_threshold_m=failure['minimum_occluded_length_m'],inset_m=1e-10),
        other_candidate_checks=dict(particle_sweep_shell_status=shell['status'],scalar_batch_status=scalar['status'],
            scalar_field_differences=scalar['differences'],max_contact_penetration_m=geometry['maximum_penetration_m']),
        interpretation=['Both centres remain in valid material space, but their straight interaction edge cuts the positive-y aperture lip.',
            'The edge already exists in the saved graphs at 30ms, 37.375ms and 39.75ms; this is not a newly introduced PID mapping error.',
            'A 0.1nm inset still leaves over 132nm of strict solid intersection, independently reproduced at 70-digit precision.',
            'Particle-path solid/contact certificates do not establish visibility of interparticle force/heat segments.',
            'The current unfiltered graph candidate is invalid under its declared visibility guard; no step from this coarse branch is accepted.'],
        limitations=['Decimal arithmetic verifies the stored endpoint geometry, not the exact continuum solution or a continuous swept graph.',
            'No individual occluded-pair force or heat effect was recomputed; a short intersection length does not justify ignoring the connection.'],
        minimum_new_model_requirements=['Apply one symmetric visibility decision to every undirected physical graph edge before any moment/correction, density, material gradient, force or thermal construction.',
            'Rebuild current geometry moments B and density-volume fixed point on the filtered graph; reject singular or ill-conditioned neighborhoods.',
            'Use the same corrected graph and B in paired internal forces and work-conjugate material gradient; preserve linear/angular momentum and virtual-work identities.',
            'Build thermal conductances and implicit graph Laplacian from the same admissible raw radial edges; retain nonnegative conductance and exact internal heat cancellation.',
            'Rebuild actual tangent/search assembly on the new graph and reevaluate the full nonlinear residual; topology/visibility changes must invalidate stale graph/tangent data.',
            'Treat visibility filtering as a new mathematical implementation with separate classification, source freeze, component identities, actual-state comparison and same-parent time refinement; never relabel this failed candidate as accepted.',
            'Declare and test the visibility event/tie/roundoff convention and endpoint-versus-swept scope; keep all existing geometric, dissipation, energy and history checks.'],
        input_sha256=inputs,source_sha256={path.relative_to(ROOT).as_posix():reader.sha(path) for path in
            (Path(__file__),Path(reader.__file__),ROOT/'scripts/audit_solid_occluded_neighbors.py',
             ROOT/'src/afsd_backend/pysph/mechanics.py',ROOT/'src/afsd_backend/pysph/quasistatic_process.py',
             ROOT/'src/afsd_backend/pysph/implicit_thermal.py',ROOT/'src/afsd_backend/pysph/quasistatic_tangent.py')})
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(status=report['status'],output=str(OUT),pids=pids,strict_nm=computed*1e9,inset_nm=inset_computed*1e9,
        roots=final['strict_intersection'],previous=[dict(time_s=item['time_s'],edge_index=item['edge_index'],
            y_at_shoulder_m=item['strict_intersection']['y_minus_positive_bore_wall_at_shoulder_intersection_m']) for item in results[:-1]])))


if __name__=='__main__':
    main()
