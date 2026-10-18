"""Conservative continuous blocked-prefix proof for declared incoming edges.

Two independently sufficient certificates keep a fixed segment parameter in
the original solid: midpoint intersection with a displacement-contracted solid,
or an endpoint-coordinate box expanded by the inherited M*h*h/8 curvature and
arithmetic bounds. Sampling alone never proves a prefix. No state, support
edge or physical gate changes.
"""
import time
import numpy as np

from visible_support import _inputs, _hash_array, _topology_hash
from visible_support_path_tight import _rotation_path
from visible_support_remote_domain_20261007 import (
    check_finite_partition_numerical_domain, require_remote_partition_consistency,
    NUMERICAL_PASS)
from audit_solid_occluded_neighbors import _partitions


CLASSIFICATION = 'CONTRACTED_SOLID_CONTINUOUS_INCOMING_BLOCKED_PREFIX_20261010_V1'
CERTIFIED = 'CERTIFIED_BLOCKED_PREFIX'
UNRESOLVED = 'UNRESOLVED'


def audit_blocked_prefix(x0, x1, i, j, pose0, pose1, p, particle_ids,
        end_fraction, maximum_queries=256, maximum_depth=32):
    """Prove every supplied segment strictly blocked for fractions [0,end].

    Particle positions and origin are linear and rotation is the original
    shortest SO(3) path. The original numerical-domain guard checks the ACTUAL
    contracted geometry before every query. A visible/degenerate sample, an
    unresolved domain, or exhausted budget returns UNRESOLVED. Invalid inputs
    raise before geometric work. This is an engineering floating-point bound,
    not a formal interval-arithmetic proof or a true-trajectory assertion.
    """
    began = time.perf_counter()
    if (type(maximum_queries) is not int or not 2 <= maximum_queries <= 256
            or type(maximum_depth) is not int or not 0 <= maximum_depth <= 32):
        raise ValueError('query budget 2..256 and subdivision depth 0..32 required')
    if (isinstance(end_fraction, (bool, np.bool_)) or not np.isscalar(end_fraction)
            or not np.isfinite(end_fraction) or not 0. < end_fraction <= 1.):
        raise ValueError('finite prefix endpoint in (0,1] required')
    end_fraction = float(end_fraction)
    first, last = np.asarray(x0, float), np.asarray(x1, float)
    if first.shape != last.shape:
        raise ValueError('matching world endpoint arrays required')
    parsed0 = _inputs(first, i, j, pose0, p, particle_ids)
    parsed1 = _inputs(last, i, j, pose1, p, particle_ids)
    _, i, j, ids, o0, R0, local0, radius, a, height, margin0 = parsed0
    _, _, _, _, o1, R1, local1, _, _, _, margin1 = parsed1
    if not len(i):
        raise ValueError('at least one incoming edge required')
    angle, rotation_at = _rotation_path(R0, R1)
    dx, do = last-first, o1-o0
    r0, r1 = first-o0, last-o1
    maxradius = np.maximum(np.linalg.norm(r0, axis=1), np.linalg.norm(r1, axis=1))
    defect = max(float(np.linalg.norm(R0.T@R0-np.eye(3), ord=np.inf)),
        float(np.linalg.norm(R1.T@R1-np.eye(3), ord=np.inf)),
        float(np.linalg.norm(rotation_at(1.)-R1, ord=np.inf)))
    scale = max(float(np.max(np.abs(first))), float(np.max(np.abs(last))),
        float(np.max(np.abs(o0))), float(np.max(np.abs(o1))), radius, a, height,
        float(maxradius.max()))
    arithmetic = 2048*np.finfo(float).eps+8*defect
    margin = float(np.nextafter(max(margin0, margin1)
        +2048*np.finfo(float).eps*scale+8*defect*float(maxradius.max()), np.inf))
    speed = (np.linalg.norm(dx-do, axis=1)+angle*maxradius)*(1.+arithmetic)+margin
    if (not np.isfinite(speed).all() or not np.isfinite(margin)
            or margin >= min(radius, a, height)/16):
        raise ArithmeticError('blocked-prefix arithmetic outside original geometric domain')
    edge_speed = np.nextafter(np.maximum(speed[i], speed[j]), np.inf)
    pairs = np.column_stack((ids[i], ids[j])).astype('<i8', copy=False)
    report = dict(status=UNRESOLVED, classification=CLASSIFICATION,
        scope='EVERY_SUPPLIED_EDGE_STRICTLY_BLOCKED_AT_EVERY_PREFIX_TIME_FRACTION',
        path='LINEAR_WORLD_PARTICLES_AND_ORIGIN_SHORTEST_SO3_TOOL_ROTATION',
        prefix_start_fraction=0., prefix_end_fraction=end_fraction,
        edge_count=len(i), maximum_queries=maximum_queries, maximum_depth=maximum_depth,
        analytic_batch_queries=0, strict_geometry_batch_queries=0, contracted_solid_queries=0,
        fixed_witness_box_queries=0, partition_queries=0,
        query_budget_scope='STRICT_PARTITION_PLUS_CONTRACTED_PARTITION_PLUS_WITNESS_BOX_BATCHES',
        deepest_subdivision=0, query_log=[], strict_samples=[], certificates=[],
        certified_edge_interval_count=0, visible_samples=[], degenerate_edges=[],
        numerical_domain_failures=[], unresolved_intervals=[], failure_kind=None,
        roundoff_margin_m=margin, rotation_angle_rad=float(angle), rotation_defect=defect,
        body_speed_bound='TRIANGLE_TRANSLATION_PLUS_ROTATION_RADIUS',
        maximum_edge_body_speed_per_fraction_m=float(edge_speed.max()),
        physical_solver_evaluations=0, nnps_evaluations=0, edge_deletions=0,
        accepted=False, restart_allowed=False,
        input_sha256=dict(world_position0=_hash_array(first), world_position1=_hash_array(last),
            particle_id_order=_hash_array(ids.astype('<i8')), incoming_edge_pid_order=_hash_array(pairs),
            incoming_topology=_topology_hash(pairs), origin0=_hash_array(o0), origin1=_hash_array(o1),
            rotation0=_hash_array(R0), rotation1=_hash_array(R1),
            geometry=_hash_array(np.array([radius,a,height])),
            prefix_endpoint=_hash_array(np.array([end_fraction]))),
        limitations=['Declared candidate chord only, not an unknown exact material trajectory.',
            'Contracted-solid intersections or conservative coordinate boxes supply fixed segment-parameter witnesses; samples alone never certify.',
            'Engineering floating-point envelopes and guarded original partitions, not formal interval arithmetic.',
            'No statement is made after the registered prefix endpoint.'])
    locals_cache = {0.:local0, 1.:local1}
    strict_cache = {}
    strict_partition_cache = {}

    class StopProof(Exception):
        pass

    def stop(reason, lo=0., hi=end_fraction, edges=None, depth=0):
        edges = np.arange(len(i)) if edges is None else edges
        report['failure_kind'] = reason
        report['unresolved_intervals'].append(dict(start_fraction=float(lo),end_fraction=float(hi),
            depth=int(depth),reason=reason,edge_indices=edges.astype(int).tolist(),
            edge_particle_ids=pairs[edges].tolist()))
        raise StopProof

    def local_at(fraction):
        if fraction not in locals_cache:
            locals_cache[fraction] = (first+fraction*dx-o0-fraction*do)@rotation_at(fraction)
        q = locals_cache[fraction]
        if not np.isfinite(q).all():
            stop('NONFINITE_BODY_PATH')
        return q

    def query(start, direction, R, A, zlo, zhi, kind, fraction, indices):
        if report['analytic_batch_queries'] >= maximum_queries:
            stop('MAXIMUM_ANALYTIC_BATCH_QUERIES')
        domain = check_finite_partition_numerical_domain(start, direction, R,
            min(radius,a,height), half_bore=A, zlo=zlo, zhi=zhi)
        if domain['status'] != NUMERICAL_PASS:
            report['numerical_domain_failures'].append(dict(kind=kind,fraction=float(fraction),
                edge_indices=indices.astype(int).tolist(),domain=domain))
            stop('PARTITION_NUMERICAL_DOMAIN_UNRESOLVED')
        report['analytic_batch_queries'] += 1
        report['partition_queries'] += 1
        report['strict_geometry_batch_queries' if kind=='strict' else 'contracted_solid_queries'] += 1
        report['query_log'].append(dict(kind=kind,fraction=float(fraction),edge_count=len(indices),
            queried_radius_m=float(R),queried_half_bore_m=float(A),queried_zlo_m=float(zlo),queried_zhi_m=float(zhi)))
        with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
            answer = _partitions(start,direction,R,A,zlo,zhi)
        if not all(np.isfinite(value).all() for value in answer):
            stop('NONFINITE_ANALYTIC_PARTITION')
        require_remote_partition_consistency(domain,answer)
        return answer

    def strict_at(fraction):
        if fraction in strict_cache:
            return
        q = local_at(fraction); direction=q[j]-q[i]
        collapsed=np.flatnonzero(np.linalg.norm(direction,axis=1)<=2*margin)
        if len(collapsed):
            report['degenerate_edges'].extend(dict(edge_index=int(edge),fraction=float(fraction),
                particle_id_i=int(ids[i[edge]]),particle_id_j=int(ids[j[edge]])) for edge in collapsed)
            stop('DEGENERATE_INCOMING_EDGE',edges=collapsed)
        cuts, inside, maximum, total = query(q[i],direction,radius,a,0.,height,'strict',fraction,np.arange(len(i)))
        report['strict_samples'].append(dict(fraction=float(fraction),
            total_inside_length_m=total.tolist(),maximum_inside_length_m=maximum.tolist()))
        strict_cache[fraction] = total.copy()
        strict_partition_cache[fraction] = (cuts,inside,maximum,total)
        visible=np.flatnonzero(total<=0.)
        if len(visible):
            report['visible_samples'].extend(dict(edge_index=int(edge),fraction=float(fraction),
                particle_id_i=int(ids[i[edge]]),particle_id_j=int(ids[j[edge]])) for edge in visible)
            stop('INCOMING_NOT_STRICTLY_BLOCKED_AT_SAMPLE',edges=visible)

    def box_witnesses(lo, hi, mid, edges):
        # This analytical proof batch shares the same <=256 overall budget;
        # it is not a free sampling/partition budget extension.
        if report['analytic_batch_queries'] >= maximum_queries:
            stop('MAXIMUM_ANALYTIC_BATCH_QUERIES')
        report['analytic_batch_queries'] += 1
        report['fixed_witness_box_queries'] += 1
        report['query_log'].append(dict(kind='fixed_witness_box',fraction=float(mid),edge_count=len(edges)))
        lower_position,upper_position=local_at(lo),local_at(hi)
        cuts,inside,_,_=strict_partition_cache[mid]
        found=[]
        for edge in edges:
            for slot in np.flatnonzero(inside[edge]):
                sa,sb=float(cuts[edge,slot]),float(cuts[edge,slot+1])
                if not sa<sb:
                    continue
                parameter=.5*(sa+sb)
                if not sa<parameter<sb:
                    continue
                n0,n1=int(i[edge]),int(j[edge])
                # The same material segment parameter is used at BOTH time
                # endpoints. Independent closest points would not be a proof.
                before=(1.-parameter)*lower_position[n0]+parameter*lower_position[n1]
                after=(1.-parameter)*upper_position[n0]+parameter*upper_position[n1]
                world0=(1.-parameter)*first[n0]+parameter*first[n1]
                world1=(1.-parameter)*last[n0]+parameter*last[n1]
                radial=max(float(np.linalg.norm(world0-o0)),float(np.linalg.norm(world1-o1)))
                velocity=float(np.linalg.norm(world1-world0-do))
                # Same SO(3) curvature bound as the frozen finite-shell box
                # audit; parameter interpolation adds the existing roundoff
                # envelope to both the derivative and coordinate enclosure.
                curvature=float(np.nextafter((2*angle*velocity+angle*angle*radial)
                    *(1.+arithmetic)+margin,np.inf))
                deviation=float(np.nextafter(curvature*(hi-lo)**2/8.+margin,np.inf))
                lower=np.nextafter(np.minimum(before,after)-deviation,-np.inf)
                upper=np.nextafter(np.maximum(before,after)+deviation,np.inf)
                radial_upper=float(np.nextafter(np.hypot(max(abs(lower[0]),abs(upper[0])),
                    max(abs(lower[1]),abs(upper[1])))*(1.+arithmetic)+margin,np.inf))
                outside_bore=any(lower[axis]>a or upper[axis]<-a for axis in (0,1))
                if (np.isfinite(lower).all() and np.isfinite(upper).all()
                        and radial_upper<radius and lower[2]>0. and upper[2]<height and outside_bore):
                    found.append(dict(edge_index=int(edge),segment_parameter_interval=[sa,sb],
                        witness_segment_parameter=parameter,body_box_lower_m=lower.tolist(),
                        body_box_upper_m=upper.tolist(),endpoint_body_positions_m=[before.tolist(),after.tolist()],
                        second_derivative_bound_m=curvature,chord_deviation_plus_roundoff_m=deviation,
                        radial_box_upper_m=radial_upper))
                    break
        return found

    try:
        # A fixed inside point is insufficient when the material segment
        # collapses between samples. World-edge vectors are affine in fraction
        # and rigid rotation preserves collapse; audit their quadratic norm
        # minimum on the complete prefix before either proof can certify it.
        direction0=first[j]-first[i]
        direction_change=(last[j]-last[i])-direction0
        for edge in range(len(i)):
            change_norm=float(np.linalg.norm(direction_change[edge]))
            if change_norm<=4*margin:
                minimum_fraction=0.
                minimum_length=float(np.linalg.norm(direction0[edge]))-change_norm*end_fraction
            else:
                unit=direction_change[edge]/change_norm
                minimum_fraction=float(np.clip(-np.dot(direction0[edge],unit)/change_norm,0.,end_fraction))
                minimum_length=float(np.linalg.norm(direction0[edge]+minimum_fraction*direction_change[edge]))
            if not np.isfinite(minimum_length) or minimum_length<=4*margin:
                report['degenerate_edges'].append(dict(edge_index=edge,fraction=minimum_fraction,
                    particle_id_i=int(ids[i[edge]]),particle_id_j=int(ids[j[edge]]),
                    minimum_world_segment_length_m=minimum_length,
                    arithmetic_exclusion_length_m=4*margin))
                stop('DEGENERATE_INCOMING_EDGE',edges=np.array([edge]))
        strict_at(0.); strict_at(end_fraction)
        stack=[(0.,end_fraction,np.arange(len(i)),0)]
        while stack:
            lo,hi,edges,depth=stack.pop()
            report['deepest_subdivision']=max(report['deepest_subdivision'],depth)
            mid=.5*(lo+hi); strict_at(mid)
            boxes=box_witnesses(lo,hi,mid,edges)
            boxed=np.array([row['edge_index'] for row in boxes],dtype=int)
            if len(boxed):
                report['certified_edge_interval_count']+=len(boxed)
                report['certificates'].append(dict(method='FIXED_SEGMENT_PARAMETER_CURVATURE_BOX',
                    start_fraction=float(lo),end_fraction=float(hi),depth=depth,reference_fraction=float(mid),
                    edge_indices=boxed.tolist(),witnesses=boxes))
            edges=np.setdiff1d(edges,boxed,assume_unique=True)
            if not len(edges):
                continue
            delta=float(np.nextafter(.5*(hi-lo)*float(edge_speed[edges].max())+margin,np.inf))
            R=float(np.nextafter(radius-delta,-np.inf))
            A=float(np.nextafter(a+delta,np.inf))
            zlo=float(np.nextafter(delta,np.inf));zhi=float(np.nextafter(height-delta,-np.inf))
            certified=np.array([],dtype=int); proof=None
            # An empty contracted solid proves nothing; retain the strict
            # sample cost and refine instead of silently clamping geometry.
            if R>0. and zlo<zhi and A<R:
                q=local_at(mid)
                proof=query(q[i[edges]],q[j[edges]]-q[i[edges]],R,A,zlo,zhi,'contracted',mid,edges)
                certified=edges[proof[3]>0.]
            if len(certified):
                witnesses=[]
                for local_index,edge in enumerate(edges):
                    if proof[3][local_index] <= 0.:
                        continue
                    slots=np.flatnonzero(proof[1][local_index])
                    slot=next(int(k) for k in slots if proof[0][local_index,k+1]>proof[0][local_index,k])
                    sa,sb=float(proof[0][local_index,slot]),float(proof[0][local_index,slot+1])
                    witnesses.append(dict(edge_index=int(edge),segment_parameter_interval=[sa,sb],
                        witness_segment_parameter=.5*(sa+sb)))
                report['certified_edge_interval_count']+=len(certified)
                report['certificates'].append(dict(method='FIXED_SEGMENT_PARAMETER_CONTRACTED_SOLID',
                    start_fraction=float(lo),end_fraction=float(hi),depth=depth,
                    reference_fraction=float(mid),displacement_plus_roundoff_bound_m=delta,
                    contracted_radius_m=R,expanded_half_bore_m=A,contracted_zlo_m=zlo,contracted_zhi_m=zhi,
                    edge_indices=certified.astype(int).tolist(),witnesses=witnesses))
            remaining=np.setdiff1d(edges,certified,assume_unique=True)
            if len(remaining):
                if depth>=maximum_depth:
                    stop('MAXIMUM_SUBDIVISION_DEPTH',lo,hi,remaining,depth)
                stack.append((mid,hi,remaining,depth+1));stack.append((lo,mid,remaining,depth+1))
        report.update(status=CERTIFIED,failure_kind=None)
    except StopProof:
        pass
    except (ValueError,ArithmeticError,FloatingPointError) as error:
        report.update(status=UNRESOLVED,failure_kind='GEOMETRIC_PROOF_ERROR',error=repr(error))
    report['wall_seconds']=time.perf_counter()-began
    return report
