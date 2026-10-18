"""Independent cancellation-aware certificate for the same fixed support graph.

The inherited triangle bound below remains a fallback. This version may use
a midpoint Rodrigues derivative plus a Frobenius-norm curvature bound; the
registered floating-point allowances enlarge the original displacement margin.
See docs/VISIBLE_SUPPORT_TIGHT_BOUND_CONTRACT_20261004.md. No physics changes.

Particles and the tool origin follow their actual declared world chords; tool
orientation follows the inherited shortest-SO(3) path. For endpoint k,
L[k] = ||dx[k]-do|| + theta*max(||r0[k]||,||r1[k]||) bounds the body-coordinate
speed per unit step fraction. Every point of edge (i,j) therefore stays within
max(L[i],L[j])*h/2 of its midpoint-time straight edge on an interval width h.

The distance-d neighborhood of the tool solid is contained in the solid with
R+d, bore half-width a-d, and z in [-d,H+d]. A reference edge disjoint from that
expanded solid proves visibility for the entire time interval, including all
edge parameters. Otherwise sample the real middle geometry and subdivide.
Sampling alone never certifies visibility; exhausted budgets are UNRESOLVED.
An interior reference edge of length <= twice the arithmetic margin is also
UNRESOLVED. The inherited length-based segment test cannot classify a collapsed
point segment, so such a degenerate graph is outside this certificate's domain.

The graph is supplied and NEVER changed by this audit. A formerly blocked but
now visible retained edge can fail at the start; this is an explicit domain
limit of the endpoint graph variant. No extra edges are deleted to obtain PASS.
"""
from pathlib import Path
import sys
import numpy as np
from visible_support import _inputs, _hash_array, _topology_hash
from audit_solid_occluded_neighbors import _partitions
from visible_support_numerics import check_partition_numerical_domain,PASS as NUMERICAL_PASS

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.quasistatic_contact import _rotation_path

CERTIFIED='CERTIFIED_VISIBLE_SUPPORT_PATH'
FOUND='SOLID_OCCLUSION_FOUND_ON_RETAINED_EDGE_PATH'
UNRESOLVED='UNRESOLVED_VISIBLE_SUPPORT_PATH'
CLASSIFICATION='FIXED_RETAINED_GRAPH_RODRIGUES_CANCELLATION_ENCLOSURE_V1'


# Independently versioned engineering enclosure; see its registered contract.
TIGHT_MINIMUM_ANGLE=1e-6
TIGHT_MAXIMUM_ANGLE=np.pi/2
TIGHT_ROUNDOFF_FACTOR=8192*np.finfo(float).eps


def _upper(value):
    return np.nextafter(np.asarray(value)*(1.+TIGHT_ROUNDOFF_FACTOR),np.inf)


def _tight_body_speed(first,last,o0,o1,R0,R1,angle,rotation_at,old_speed,old_margin,scale):
    """Bound the *same* frozen Rodrigues polynomial, preserving cancellation.

    The exact reference polynomial uses the stored R0, K, K@K coefficients;
    it does not assume they satisfy exact orthogonality identities. The
    engineering floating-point allowances cover evaluation and cancellation.
    No geometric queries, edges, coordinates, or tool poses are changed.
    """
    info=dict(method='INHERITED_TRIANGLE_BOUND_FALLBACK',minimum_tight_angle_rad=TIGHT_MINIMUM_ANGLE,
        maximum_tight_angle_rad=float(TIGHT_MAXIMUM_ANGLE),roundoff_factor=TIGHT_ROUNDOFF_FACTOR,
        original_displacement_roundoff_margin_m=float(old_margin),extra_displacement_roundoff_margin_m=0.,
        maximum_original_body_speed_per_fraction_m=float(np.max(old_speed,initial=0.)),
        tightened_particle_count=0,fallback_reason=None)
    if not TIGHT_MINIMUM_ANGLE<=angle<=TIGHT_MAXIMUM_ANGLE:
        info['fallback_reason']='ROTATION_ANGLE_OUTSIDE_REGISTERED_TIGHT_DOMAIN'
        return old_speed.copy(),old_margin,info
    # These operations match the frozen _rotation_path axis extraction.
    relative=R0.T@R1
    skew=.5*np.array([relative[2,1]-relative[1,2],relative[0,2]-relative[2,0],relative[1,0]-relative[0,1]])
    sine=float(np.linalg.norm(skew))
    if not np.isfinite(sine) or sine<=0.:raise ArithmeticError('unresolved tight-path rotation axis')
    axis=skew/sine
    K=np.array([[0.,-axis[2],axis[1]],[axis[2],0.,-axis[0]],[-axis[1],axis[0],0.]])
    K2=K@K
    nR=float(_upper(np.linalg.norm(R0,ord='fro')))
    nK=float(_upper(np.linalg.norm(K,ord='fro')))
    nK2=float(_upper(np.linalg.norm(K2,ord='fro')))
    C=float(_upper(nR*(nK+nK2)))
    A=float(_upper(nR*(float(_upper(np.sqrt(3.)))+nK+2*nK2)))
    dx=last-first;do=o1-o0;v=dx-do
    rho=_upper(np.maximum(np.linalg.norm(first-o0,axis=1),np.linalg.norm(last-o1,axis=1)))
    vn=_upper(np.linalg.norm(v,axis=1))
    mid_r=first+.5*dx-o0-.5*do
    Rmid=rotation_at(.5)
    Rprime=R0@(angle*(np.cos(.5*angle)*K+np.sin(.5*angle)*K2))
    translation_term=v@Rmid;rotation_term=mid_r@Rprime
    derivative=np.linalg.norm(translation_term+rotation_term,axis=1)
    curvature=_upper(C*(2*angle*vn+angle*angle*rho))
    derivative_error=_upper(TIGHT_ROUNDOFF_FACTOR*(A*(scale+vn)+angle*C*(scale+rho)))
    tight=_upper(derivative+.5*curvature+derivative_error)
    # Uniform evaluation error at both moving point and midpoint reference.
    path_error=float(_upper(TIGHT_ROUNDOFF_FACTOR*A*(scale+float(rho.max(initial=0.)))))
    margin=float(np.nextafter(old_margin+2*path_error,np.inf))
    chosen=np.minimum(old_speed,tight)
    if not (np.isfinite(chosen).all() and np.isfinite(tight).all() and
            np.isfinite(curvature).all() and np.isfinite(derivative_error).all() and np.isfinite(margin)):
        raise ArithmeticError('nonfinite tight-path bound')
    info.update(method='MIN_INHERITED_AND_MIDPOINT_DERIVATIVE_PLUS_HALF_CURVATURE',
        extra_displacement_roundoff_margin_m=float(margin-old_margin),
        uniform_position_evaluation_error_m=path_error,
        rodrigues_derivative_coefficient_upper=C,rodrigues_matrix_norm_upper=A,
        tightened_particle_count=int(np.count_nonzero(tight<old_speed)),
        maximum_tight_body_speed_per_fraction_m=float(tight.max(initial=0.)),
        maximum_selected_body_speed_per_fraction_m=float(chosen.max(initial=0.)),
        per_particle_old_speed_m=old_speed.tolist(),per_particle_tight_speed_m=tight.tolist(),
        per_particle_selected_speed_m=chosen.tolist(),per_particle_midpoint_derivative_m=derivative.tolist(),
        per_particle_curvature_upper_m=curvature.tolist(),per_particle_absolute_derivative_roundoff_m=derivative_error.tolist(),
        limitations=['Engineering binary64 rounding envelope, not formal interval arithmetic.',
            'The small midpoint derivative is not assigned a relative-error-only allowance.'])
    return chosen,margin,info


def audit_visible_support_path(x0,x1,i,j,pose0,pose1,p,particle_ids,maximum_depth=32,maximum_queries=256):
    """Certify only the given fixed edges on the declared whole-interval path.

    A query is one analytic batch segment/solid intersection call, including
    strict endpoint/midpoint calls AND inflated-solid proof calls. Both count
    against the same <=256 budget. No force, NNPS, thermal or material work.
    """
    if (isinstance(maximum_depth,bool) or not isinstance(maximum_depth,int) or not 0<=maximum_depth<=32 or
        isinstance(maximum_queries,bool) or not isinstance(maximum_queries,int) or not 2<=maximum_queries<=256):
        raise ValueError('depth 0..32 and analytic batch query budget 2..256 required')
    first=np.asarray(x0,float);last=np.asarray(x1,float)
    if first.shape!=last.shape:raise ValueError('matching world endpoint position arrays required')
    parsed0=_inputs(first,i,j,pose0,p,particle_ids);parsed1=_inputs(last,i,j,pose1,p,particle_ids)
    _,i,j,ids,o0,R0,local0,radius,a,height,margin0=parsed0
    _,_,_,_,o1,R1,local1,_,_,_,margin1=parsed1
    angle,rotation_at=_rotation_path(R0,R1)
    dx=last-first;do=o1-o0;r0=first-o0;r1=last-o1
    maxradius=np.maximum(np.linalg.norm(r0,axis=1),np.linalg.norm(r1,axis=1))
    defect=max(float(np.linalg.norm(R0.T@R0-np.eye(3),ord=np.inf)),
        float(np.linalg.norm(R1.T@R1-np.eye(3),ord=np.inf)),float(np.linalg.norm(rotation_at(1.)-R1,ord=np.inf)))
    scale=max(float(np.max(np.abs(first))),float(np.max(np.abs(last))),
        float(np.max(np.abs(o0))),float(np.max(np.abs(o1))),radius,a,height,float(maxradius.max(initial=0.)))
    margin=float(np.nextafter(max(margin0,margin1)+2048*np.finfo(float).eps*scale+8*defect*float(maxradius.max(initial=0.)),np.inf))
    arithmetic=2048*np.finfo(float).eps+8*defect
    speed=(np.linalg.norm(dx-do,axis=1)+angle*maxradius)*(1.+arithmetic)+margin
    speed,margin,speed_bound_report=_tight_body_speed(first,last,o0,o1,R0,R1,angle,rotation_at,speed,margin,scale)
    if margin>=min(radius,a,height)/16:raise ArithmeticError('tight-path arithmetic margin outside registered geometric domain')
    edge_speed=np.maximum(speed[i],speed[j])
    pid_pairs=np.column_stack((ids[i],ids[j])).astype('<i8',copy=False)
    report=dict(status=UNRESOLVED,classification=CLASSIFICATION,
        scope='SUPPLIED_FIXED_RETAINED_ENDPOINT_GRAPH_ALL_TIME_FRACTIONS_AND_ALL_EDGE_PARAMETERS',
        path='LINEAR_WORLD_PARTICLES_AND_ORIGIN_SHORTEST_SO3_TOOL_ROTATION',
        physical_solver_evaluations=0,nnps_evaluations=0,edge_deletions=0,edge_count=len(i),
        maximum_queries=maximum_queries,maximum_depth=maximum_depth,analytic_batch_queries=0,
        strict_geometry_batch_queries=0,enclosure_batch_queries=0,deepest_subdivision=0,
        evaluated_fractions=[],query_log=[],certified_edge_interval_count=0,certificates=[],
        found_edges=[],degenerate_edges=[],numerical_domain_failures=[],unresolved_intervals=[],failure_kind=None,
        rotation_angle_rad=float(angle),rotation_defect=defect,roundoff_margin_m=margin,
        speed_bound=speed_bound_report,
        maximum_edge_body_speed_per_fraction_m=float(edge_speed.max(initial=0.)),
        input_sha256=dict(world_position0=_hash_array(first),world_position1=_hash_array(last),
            particle_id_order=_hash_array(ids.astype('<i8')),retained_edge_pid_order=_hash_array(pid_pairs),
            retained_topology=_topology_hash(pid_pairs),origin0=_hash_array(o0),origin1=_hash_array(o1),
            rotation0=_hash_array(R0),rotation1=_hash_array(R1),
            geometry=_hash_array(np.array([radius,a,height]))),
        limitations=['Only the supplied declared path is certified, not an unknown exact material trajectory.',
            'A start-blocked/end-visible kept edge causes rejection; this audit never deletes more edges.',
            'Collapsed or arithmetic-scale interior support edges are unresolved, even if the collapsed point is outside the tool.',
            'Floating-point arithmetic with explicit conservative margins, not a formal interval-arithmetic proof.'])
    locals_cache={0.:local0,1.:local1};strict_cache={}

    class Budget(Exception):pass
    class Degenerate(Exception):pass
    class NumericalDomain(Exception):pass

    def query(start,direction,R,A,zlo,zhi,kind,fraction,count):
        if report['analytic_batch_queries']>=maximum_queries:raise Budget
        domain=check_partition_numerical_domain(start,direction,R,min(radius,a,height))
        if domain['status']!=NUMERICAL_PASS:
            report['numerical_domain_failures'].append(dict(kind=kind,fraction=float(fraction),edge_count=int(count),domain=domain))
            unresolved(0.,1.,np.arange(len(i)),0,'PARTITION_NUMERICAL_DOMAIN_UNRESOLVED')
            raise NumericalDomain
        report['analytic_batch_queries']+=1
        report['strict_geometry_batch_queries' if kind=='strict' else 'enclosure_batch_queries']+=1
        report['query_log'].append(dict(kind=kind,fraction=float(fraction),edge_count=int(count)))
        with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
            result=_partitions(start,direction,R,A,zlo,zhi)
        if not all(np.isfinite(value).all() for value in result):raise ArithmeticError('nonfinite segment/path enclosure')
        return result

    def local_at(fraction):
        if fraction not in locals_cache:
            rotation=rotation_at(fraction)
            locals_cache[fraction]=(first+fraction*dx-o0-fraction*do)@rotation
        return locals_cache[fraction]

    def strict_at(fraction):
        if fraction not in strict_cache:
            q=local_at(fraction)
            direction=q[j]-q[i]
            collapsed=np.flatnonzero(np.linalg.norm(direction,axis=1)<=2*margin)
            if len(collapsed):
                for edge in collapsed:
                    report['degenerate_edges'].append(dict(edge_index=int(edge),i=int(i[edge]),j=int(j[edge]),
                        particle_id_i=int(ids[i[edge]]),particle_id_j=int(ids[j[edge]]),fraction=float(fraction),
                        reference_edge_length_m=float(np.linalg.norm(direction[edge])),
                        body_start=q[i[edge]].tolist(),body_end=q[j[edge]].tolist()))
                unresolved(0.,1.,collapsed,0,'DEGENERATE_INTERIOR_SUPPORT_EDGE')
                raise Degenerate
            cuts,inside,maximum,total=query(q[i],direction,radius,a,0.,height,'strict',fraction,len(i))
            report['evaluated_fractions'].append(float(fraction))
            strict_cache[fraction]=total
            for edge in np.flatnonzero(total>0.):
                report['found_edges'].append(dict(edge_index=int(edge),i=int(i[edge]),j=int(j[edge]),
                    particle_id_i=int(ids[i[edge]]),particle_id_j=int(ids[j[edge]]),fraction=float(fraction),
                    maximum_contiguous_inside_length_m=float(maximum[edge]),total_inside_length_m=float(total[edge]),
                    strict_intervals=[[float(cuts[edge,n]),float(cuts[edge,n+1])] for n in range(inside.shape[1])
                        if inside[edge,n] and cuts[edge,n+1]>cuts[edge,n]]))
        return strict_cache[fraction]

    def unresolved(lo,hi,edges,depth,reason):
        report['unresolved_intervals'].append(dict(start_fraction=float(lo),end_fraction=float(hi),depth=int(depth),
            reason=reason,edge_indices=edges.astype(int).tolist(),edge_particle_ids=pid_pairs[edges].tolist()))

    stack=[]
    try:
        strict_at(0.);strict_at(1.)
        if not report['found_edges'] and len(i):stack=[(0.,1.,np.arange(len(i)),0)]
        while stack and not report['found_edges']:
            lo,hi,edges,depth=stack.pop();report['deepest_subdivision']=max(report['deepest_subdivision'],depth)
            mid=.5*(lo+hi)
            try:
                strict_at(mid)
                if report['found_edges']:break
                q=local_at(mid)
                displacement_bound=float(np.nextafter(.5*(hi-lo)*float(edge_speed[edges].max())+margin,np.inf))
                # a-d may be negative: the expanded solid then has no hole.
                # Do not clamp it to zero, which would accidentally spare the
                # reference segment lying exactly on the tool axis.
                _,_,_,total=query(q[i[edges]],q[j[edges]]-q[i[edges]],
                    np.nextafter(radius+displacement_bound,np.inf),np.nextafter(a-displacement_bound,-np.inf),
                    np.nextafter(-displacement_bound,-np.inf),np.nextafter(height+displacement_bound,np.inf),
                    'enclosure',mid,len(edges))
            except Budget:
                unresolved(lo,hi,edges,depth,'MAXIMUM_ANALYTIC_BATCH_QUERIES')
                for left,right,remaining,level in stack:unresolved(left,right,remaining,level,'MAXIMUM_ANALYTIC_BATCH_QUERIES')
                stack=[];report['failure_kind']='MAXIMUM_ANALYTIC_BATCH_QUERIES';break
            certified=edges[total==0.];remaining=edges[total>0.]
            if len(certified):
                report['certified_edge_interval_count']+=len(certified)
                report['certificates'].append(dict(start_fraction=float(lo),end_fraction=float(hi),
                    depth=depth,reference_fraction=float(mid),displacement_plus_roundoff_bound_m=displacement_bound,
                    edge_indices=certified.astype(int).tolist()))
            if len(remaining):
                if depth>=maximum_depth:
                    unresolved(lo,hi,remaining,depth,'MAXIMUM_SUBDIVISION_DEPTH')
                    report['failure_kind']='MAXIMUM_SUBDIVISION_DEPTH'
                else:
                    stack.append((mid,hi,remaining,depth+1));stack.append((lo,mid,remaining,depth+1))
    except Budget:
        unresolved(0.,1.,np.arange(len(i)),0,'MAXIMUM_ANALYTIC_BATCH_QUERIES')
        report['failure_kind']='MAXIMUM_ANALYTIC_BATCH_QUERIES'
    except Degenerate:
        report['failure_kind']='DEGENERATE_INTERIOR_SUPPORT_EDGE'
    except NumericalDomain:
        report['failure_kind']='PARTITION_NUMERICAL_DOMAIN_UNRESOLVED'
    if report['found_edges']:
        report['status']=FOUND;report['failure_kind']='STRICT_SOLID_INTERSECTION_ON_FIXED_RETAINED_PATH'
    elif not report['unresolved_intervals']:
        report['status']=CERTIFIED;report['failure_kind']=None
    return report
