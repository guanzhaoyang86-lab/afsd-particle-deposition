"""Conservative continuous visibility using inherited proof plus curvature prisms.

Every original strict/domain and inflated-solid query remains in place. Only
unresolved edges can additionally be proved disjoint from four wall-rectangle
supersets by a four-endpoint convex hull plus a curvature/roundoff box. Each
SAT batch shares the original 256-query limit and 32 subdivision depth.
These are engineering binary64 enclosures, not formal interval arithmetic.
No support edge, physical state, time interval or acceptance gate is changed.
"""
import numpy as np
from visible_support import _inputs, _hash_array, _topology_hash
from audit_solid_occluded_neighbors import _partitions
from visible_support_remote_domain_20261007 import (
    check_finite_partition_numerical_domain, require_remote_partition_consistency,
    NUMERICAL_PASS, CERTIFIED, FOUND, UNRESOLVED, _rotation_path, _tight_body_speed)

CLASSIFICATION='FIXED_RETAINED_GRAPH_CURVATURE_PRISM_SAT_20261010_V1'
FROZEN_PARENT_SHA256='197a2f57fc9b33712e5858e08d036df9b0c2bbfcda8814d1d8781aa74c310067'
PROJECTION_ROUNDOFF_FACTOR=2048*np.finfo(float).eps


def _curvature_bounds(first,last,o0,o1,R0,R1,angle,speed_report):
    """Second derivatives use the original complete path fraction variable."""
    if 'per_particle_curvature_upper_m' in speed_report:
        curvature=np.asarray(speed_report['per_particle_curvature_upper_m'],float)
        method='INHERITED_RODRIGUES_FROBENIUS_SECOND_DERIVATIVE'
    elif angle==0.:
        curvature=np.zeros(len(first));method='CONSTANT_ROTATION_AFFINE_BODY_PATH'
    else:
        return None,'SAT_UNAVAILABLE_OUTSIDE_INHERITED_CURVATURE_DOMAIN'
    if curvature.shape!=(len(first),) or not np.isfinite(curvature).all() or np.any(curvature<0):
        raise ArithmeticError('nonfinite/negative prism curvature bound')
    return curvature,method


def _separate_four_wall_rectangles(points,delta,radius,a,height,margin):
    """All four signed-horizontal/z rectangle supersets must strictly separate.

    Every finite nonzero projection axis is a valid separation test; the six
    pair normals merely propose useful axes and need not be exact hull normals.
    The unknown other horizontal coordinate is never constrained.
    """
    points=np.asarray(points,float)
    if points.shape!=(4,3) or not np.isfinite(points).all() or not np.isfinite(delta) or delta<margin:
        return None
    walls=[]
    rect=np.array([[a,0.],[a,height],[radius,0.],[radius,height]])
    for axis in (0,1):
        for sign in (-1,1):
            planar=np.column_stack((sign*points[:,axis],points[:,2]))
            axes=[np.array([1.,0.]),np.array([0.,1.])]
            for left in range(4):
                for right in range(left+1,4):
                    d=planar[right]-planar[left]
                    magnitude=float(np.max(np.abs(d)))
                    if np.isfinite(magnitude) and magnitude>4*margin:
                        axes.append(np.array([-d[1],d[0]])/magnitude)
            witness=None
            for normal in axes:
                if not np.isfinite(normal).all() or np.max(np.abs(normal))==0.:continue
                l1=float(np.nextafter(np.sum(np.abs(normal)),np.inf))
                point_projection=planar@normal;rectangle_projection=rect@normal
                # Each dot-product support interval is enlarged separately,
                # before testing the gap. Cancellation never sets the scale.
                point_scale=float(np.max(np.sum(np.abs(planar)*np.abs(normal),axis=1)))
                rectangle_scale=float(np.max(np.sum(np.abs(rect)*np.abs(normal),axis=1)))
                point_error=float(np.nextafter(PROJECTION_ROUNDOFF_FACTOR*point_scale,np.inf))
                rectangle_error=float(np.nextafter(PROJECTION_ROUNDOFF_FACTOR*rectangle_scale,np.inf))
                padding=float(np.nextafter(np.nextafter(delta*l1,np.inf)+point_error,np.inf))
                hull_lo=float(np.nextafter(point_projection.min()-padding,-np.inf))
                hull_hi=float(np.nextafter(point_projection.max()+padding,np.inf))
                rect_lo=float(np.nextafter(rectangle_projection.min()-rectangle_error,-np.inf))
                rect_hi=float(np.nextafter(rectangle_projection.max()+rectangle_error,np.inf))
                if not np.isfinite([hull_lo,hull_hi,rect_lo,rect_hi]).all():continue
                if hull_hi<rect_lo:
                    orientation='HULL_BELOW_RECTANGLE';gap=float(np.nextafter(rect_lo-hull_hi,-np.inf))
                elif rect_hi<hull_lo:
                    orientation='RECTANGLE_BELOW_HULL';gap=float(np.nextafter(hull_lo-rect_hi,-np.inf))
                else:continue
                if gap<=0.:continue
                witness=dict(horizontal_axis='xy'[axis],horizontal_sign=sign,normal=normal.tolist(),
                    hull_projection_interval_m=[hull_lo,hull_hi],rectangle_projection_interval_m=[rect_lo,rect_hi],
                    curvature_roundoff_projection_padding_m=padding,point_dot_roundoff_m=point_error,
                    rectangle_dot_roundoff_m=rectangle_error,strict_gap_m=gap,orientation=orientation)
                break
            if witness is None:return None
            walls.append(witness)
    return walls


def audit_visible_support_path(x0,x1,i,j,pose0,pose1,p,particle_ids,maximum_depth=32,maximum_queries=256):
    """Certify only the given fixed edges on the declared whole-interval path.

    A query is one strict partition, inflated-solid partition, or additional
    prism SAT batch. All count against the same <=256 budget. No force, NNPS,
    thermal or material work is performed.
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
    curvature,curvature_method=_curvature_bounds(first,last,o0,o1,R0,R1,angle,speed_bound_report)
    pid_pairs=np.column_stack((ids[i],ids[j])).astype('<i8',copy=False)
    report=dict(status=UNRESOLVED,classification=CLASSIFICATION,
        scope='SUPPLIED_FIXED_RETAINED_ENDPOINT_GRAPH_ALL_TIME_FRACTIONS_AND_ALL_EDGE_PARAMETERS',
        path='LINEAR_WORLD_PARTICLES_AND_ORIGIN_SHORTEST_SO3_TOOL_ROTATION',
        physical_solver_evaluations=0,nnps_evaluations=0,edge_deletions=0,edge_count=len(i),
        maximum_queries=maximum_queries,maximum_depth=maximum_depth,analytic_batch_queries=0,
        strict_geometry_batch_queries=0,enclosure_batch_queries=0,prism_batch_queries=0,deepest_subdivision=0,
        frozen_parent_sha256=FROZEN_PARENT_SHA256,
        query_budget_scope='STRICT_PARTITION_PLUS_INFLATED_PARTITION_PLUS_PRISM_SAT_BATCHES',
        evaluated_fractions=[],query_log=[],certified_edge_interval_count=0,certificates=[],
        found_edges=[],degenerate_edges=[],numerical_domain_failures=[],unresolved_intervals=[],failure_kind=None,
        rotation_angle_rad=float(angle),rotation_defect=defect,roundoff_margin_m=margin,
        speed_bound=speed_bound_report,prism_curvature_method=curvature_method,
        prism_scope='ALL_EDGE_PARAMETERS_IN_FOUR_ENDPOINT_CONVEX_HULL_PLUS_CURVATURE_BOX',
        maximum_edge_body_speed_per_fraction_m=float(edge_speed.max(initial=0.)),
        input_sha256=dict(world_position0=_hash_array(first),world_position1=_hash_array(last),
            particle_id_order=_hash_array(ids.astype('<i8')),retained_edge_pid_order=_hash_array(pid_pairs),
            retained_topology=_topology_hash(pid_pairs),origin0=_hash_array(o0),origin1=_hash_array(o1),
            rotation0=_hash_array(R0),rotation1=_hash_array(R1),
            geometry=_hash_array(np.array([radius,a,height]))),
        limitations=['Only the supplied declared path is certified, not an unknown exact material trajectory.',
            'A start-blocked/end-visible kept edge causes rejection; this audit never deletes more edges.',
            'Collapsed or arithmetic-scale interior support edges are unresolved, even if the collapsed point is outside the tool.',
            'Floating-point arithmetic with explicit conservative margins, not a formal interval-arithmetic proof.',
            'Each of four signed-horizontal/z wall rectangles is a superset; the other horizontal coordinate is unconstrained.'])
    locals_cache={0.:local0,1.:local1};strict_cache={}

    class Budget(Exception):pass
    class Degenerate(Exception):pass
    class NumericalDomain(Exception):pass

    def query(start,direction,R,A,zlo,zhi,kind,fraction,count):
        if report['analytic_batch_queries']>=maximum_queries:raise Budget
        domain=check_finite_partition_numerical_domain(start,direction,R,min(radius,a,height), half_bore=A, zlo=zlo, zhi=zhi)
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
        require_remote_partition_consistency(domain, result)
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

    def prism_witnesses(lo,hi,edges):
        if curvature is None:return []
        if report['analytic_batch_queries']>=maximum_queries:raise Budget
        report['analytic_batch_queries']+=1;report['prism_batch_queries']+=1
        report['query_log'].append(dict(kind='prism_sat',fraction=float(.5*(lo+hi)),edge_count=len(edges)))
        lower,upper=local_at(lo),local_at(hi)
        found=[]
        for edge in edges:
            points=np.array([lower[i[edge]],lower[j[edge]],upper[i[edge]],upper[j[edge]]])
            second=float(max(curvature[i[edge]],curvature[j[edge]]))
            delta=float(np.nextafter(np.nextafter(second*(hi-lo)**2/8.,np.inf)+margin,np.inf))
            walls=_separate_four_wall_rectangles(points,delta,radius,a,height,margin)
            if walls is not None:
                found.append(dict(edge_index=int(edge),body_endpoint_hull_points_m=points.tolist(),
                    second_derivative_bound_m=second,curvature_plus_roundoff_m=delta,
                    wall_separations=walls))
        return found

    stack=[]
    try:
        # A four-point hull can enclose a collapsed edge outside the solid.
        # This additional guard enforces the nondegenerate-edge domain over
        # ALL fractions, including unsampled times.
        direction0=first[j]-first[i]
        change=(last[j]-last[i])-direction0
        change_norm=np.linalg.norm(change,axis=1)
        moving=change_norm>4*margin
        fraction=np.zeros(len(i))
        fraction[moving]=np.clip(-np.einsum('ij,ij->i',direction0[moving],
            change[moving]/change_norm[moving,None])/change_norm[moving],0.,1.)
        lengths=np.linalg.norm(direction0+fraction[:,None]*change,axis=1)
        lengths[~moving]=np.linalg.norm(direction0[~moving],axis=1)-change_norm[~moving]
        collapsed=np.flatnonzero(~np.isfinite(lengths)|(lengths<=4*margin))
        if len(collapsed):
            for edge in collapsed:
                report['degenerate_edges'].append(dict(edge_index=int(edge),i=int(i[edge]),j=int(j[edge]),
                    particle_id_i=int(ids[i[edge]]),particle_id_j=int(ids[j[edge]]),fraction=float(fraction[edge]),
                    reference_edge_length_m=float(lengths[edge]),arithmetic_exclusion_length_m=4*margin))
            unresolved(0.,1.,collapsed,0,'DEGENERATE_INTERIOR_SUPPORT_EDGE')
            raise Degenerate
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
                prisms=prism_witnesses(lo,hi,edges[total>0.]) if np.any(total>0.) else []
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
            if prisms:
                prism_edges=np.array([row['edge_index'] for row in prisms],dtype=int)
                report['certified_edge_interval_count']+=len(prism_edges)
                report['certificates'].append(dict(method='FOUR_WALL_RECTANGLE_CURVATURE_PRISM_SAT',
                    start_fraction=float(lo),end_fraction=float(hi),depth=depth,
                    edge_indices=prism_edges.tolist(),witnesses=prisms))
                remaining=np.setdiff1d(remaining,prism_edges,assume_unique=True)
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
