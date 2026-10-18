"""Finite-tool broad phase before the unchanged near-tool numerical guard.

Only segments proved outside the actual queried finite cylinder are exempted
from the inherited cylinder-root numerical-domain check.  The same full edge
arrays still reach the original filter, force and thermal operators.  Explicit
postconditions reject any disagreement between the proof and that filter or
partitioner; no support edge is silently dropped or added.

The outward bounds below use declared binary64 engineering roundoff envelopes,
not formal interval arithmetic.  Near and ambiguous segments retain the exact
old domain check, including its eight-tool-scale and discriminant tests.
"""
import ast
import inspect
import textwrap
import numpy as np

from visible_support import _inputs, _hash_array, _topology_hash
from visible_support_numerics import (
    check_partition_numerical_domain as legacy_check,
    VisibleSupportNumericalDomainError, PASS, UNRESOLVED as DOMAIN_UNRESOLVED,
)
from audit_solid_occluded_neighbors import _partitions
import visible_support_path_tight as legacy_path

DOMAIN_CLASSIFICATION = 'FINITE_TOOL_PROVED_REMOTE_DOMAIN_ROUTING_20261007'
CLASSIFICATION = 'FIXED_RETAINED_GRAPH_TIGHT_PATH_FINITE_TOOL_REMOTE_DOMAIN_V1'
CERTIFIED, FOUND, UNRESOLVED = legacy_path.CERTIFIED, legacy_path.FOUND, legacy_path.UNRESOLVED
NUMERICAL_PASS = PASS
_rotation_path = legacy_path._rotation_path
_tight_body_speed = legacy_path._tight_body_speed
BROAD_PHASE_ROUNDOFF_FACTOR = 64*np.finfo(float).eps
# The unchanged partitioner forms fourth-degree discriminant products.  Remote
# geometry does not license overflow in a call that still receives those rows.
MAXIMUM_PARTITION_MAGNITUDE = np.finfo(float).max**.25/16.


def check_finite_partition_numerical_domain(start, direction, radius,
        minimum_tool_dimension, *, half_bore, zlo, zhi):
    """Audit one exact strict or inflated partition query without changing rows.

    The supplied R, a, zlo and zhi must be the ACTUAL query dimensions.  A
    negative half_bore is legitimate for an inflated solid with no open hole.
    The bore is not needed by the outside-cylinder proof.  The distance from
    the axis to the segment's outward XY bounding box is a conservative lower
    bound for the distance to the whole segment; endpoints alone are not used.
    """
    s, d = np.asarray(start, float), np.asarray(direction, float)
    radius, scale = float(radius), float(minimum_tool_dimension)
    half_bore, zlo, zhi = float(half_bore), float(zlo), float(zhi)
    if (s.ndim != 2 or s.shape[1] != 3 or d.shape != s.shape
            or not np.isfinite(s).all() or not np.isfinite(d).all()):
        raise ValueError('finite matching N by 3 segment arrays required')
    if (not np.isfinite([radius, scale, half_bore, zlo, zhi]).all()
            or min(radius, scale) <= 0 or zlo >= zhi):
        raise ValueError('finite positive queried radius/original scale and ordered slab required')
    n = len(s)
    report = dict(status=PASS, classification=DOMAIN_CLASSIFICATION,
        scope='PROVED_OUTSIDE_ACTUAL_QUERIED_FINITE_CYLINDER_OR_UNCHANGED_LEGACY_GUARD',
        edge_count=n, queried_cylinder_radius_m=radius, queried_half_bore_m=half_bore,
        queried_zlo_m=zlo, queried_zhi_m=zhi, original_minimum_tool_dimension_m=scale,
        reasons=[], uncertain_edge_indices=[], remote_edge_indices=[],
        remote_z_edge_indices=[], remote_radial_edge_indices=[], near_edge_indices=[],
        broad_phase_roundoff_factor=BROAD_PHASE_ROUNDOFF_FACTOR,
        maximum_partition_magnitude=MAXIMUM_PARTITION_MAGNITUDE,
        edge_deletions=0, graph_updates=0, force_evaluations=0,
        physical_solver_evaluations=0, nnps_evaluations=0,
        legacy_domain=None, proof='OUTWARD_SEGMENT_AABB_OUTSIDE_QUERY_Z_SLAB_OR_CYLINDER',
        limitations=['Declared engineering floating-point envelope, not formal interval arithmetic.',
            'A query-local proof does not certify a moving segment over a time interval.',
            'Every proved remote row remains in the original partition and SPH operators.'])
    magnitude = max(float(np.max(np.abs(s), initial=0.)),
        float(np.max(np.abs(d), initial=0.)), radius, abs(half_bore), abs(zlo), abs(zhi))
    if magnitude > MAXIMUM_PARTITION_MAGNITUDE:
        report.update(status=DOMAIN_UNRESOLVED,
            reasons=['UNCHANGED_PARTITION_ARITHMETIC_MAGNITUDE_UNRESOLVED'],
            uncertain_edge_indices=list(range(n)))
        return report
    with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
        end = s+d
        # Covers endpoint addition and the original partitioner's intermediate
        # affine sample evaluations, including subtraction at a slab boundary.
        padding = np.nextafter(BROAD_PHASE_ROUNDOFF_FACTOR*(np.abs(s)+np.abs(d)
            +max(radius, abs(zlo), abs(zhi))), np.inf)
        lower = np.nextafter(np.minimum(s, np.nextafter(end, -np.inf))-padding, -np.inf)
        upper = np.nextafter(np.maximum(s, np.nextafter(end, np.inf))+padding, np.inf)
        remote_z = (upper[:, 2] < zlo) | (lower[:, 2] > zhi)
        axis_gap = np.where(lower[:, :2] > 0., lower[:, :2],
            np.where(upper[:, :2] < 0., -upper[:, :2], 0.))
        squared_lower = np.maximum(0., np.nextafter(axis_gap*axis_gap, -np.inf))
        distance_squared_lower = np.maximum(0., np.nextafter(
            squared_lower[:, 0]+squared_lower[:, 1], -np.inf))
        radius_squared_upper = np.nextafter(radius*radius, np.inf)
        remote_radial = distance_squared_lower > radius_squared_upper
    remote = remote_z | remote_radial
    near_indices = np.flatnonzero(~remote)
    report.update(remote_edge_indices=np.flatnonzero(remote).astype(int).tolist(),
        remote_z_edge_indices=np.flatnonzero(remote_z).astype(int).tolist(),
        remote_radial_edge_indices=np.flatnonzero(remote_radial).astype(int).tolist(),
        near_edge_indices=near_indices.astype(int).tolist())
    if len(near_indices):
        old = legacy_check(s[near_indices], d[near_indices], radius, scale)
        report['legacy_domain'] = old
        report['status'] = old['status']
        report['reasons'] = list(old['reasons'])
        report['uncertain_edge_indices'] = [int(near_indices[k])
            for k in old['uncertain_edge_indices']]
    return report


def require_remote_partition_consistency(domain, result):
    """Check the unchanged partition result; disagreement rejects the query."""
    if domain.get('classification') != DOMAIN_CLASSIFICATION or domain['status'] != PASS:
        raise ValueError('successful finite-tool domain report required')
    remote = np.asarray(domain['remote_edge_indices'], dtype=int)
    cuts, inside, maximum, total = result
    if (np.asarray(cuts).shape[0] != domain['edge_count']
            or np.asarray(inside).shape[0] != domain['edge_count']
            or np.asarray(maximum).shape != (domain['edge_count'],)
            or np.asarray(total).shape != (domain['edge_count'],)):
        raise ValueError('partition result does not match audited rows')
    if (np.any(np.asarray(inside)[remote]) or np.any(np.asarray(maximum)[remote] != 0.)
            or np.any(np.asarray(total)[remote] != 0.)):
        raise ArithmeticError('unchanged partition disagrees with proved remote support')


def require_visible_support_numerical_domain(x, i, j, pose, p, particle_ids):
    """Same endpoint API, proving strict AND actual expanded query separation."""
    _, i, j, _, _, _, local, radius, a, height, margin = _inputs(x, i, j, pose, p, particle_ids)
    start, direction = local[i], local[j]-local[i]
    scale = min(radius, a, height)
    strict = check_finite_partition_numerical_domain(start, direction, radius, scale,
        half_bore=a, zlo=0., zhi=height)
    expanded = check_finite_partition_numerical_domain(start, direction,
        np.nextafter(radius+margin, np.inf), scale,
        half_bore=np.nextafter(a-margin, -np.inf), zlo=np.nextafter(-margin, -np.inf),
        zhi=np.nextafter(height+margin, np.inf))
    report = dict(status=PASS if strict['status'] == PASS and expanded['status'] == PASS
        else DOMAIN_UNRESOLVED, classification=DOMAIN_CLASSIFICATION,
        scope='STRICT_AND_ACTUAL_ROUNDOFF_EXPANDED_FINITE_TOOL_ENDPOINT_QUERIES',
        strict=strict, expanded=expanded,
        reasons=sorted(set(strict['reasons']+expanded['reasons'])),
        proven_visible_raw_edge_indices=sorted(set(strict['remote_edge_indices'])
            & set(expanded['remote_edge_indices'])),
        edge_count=len(i), edge_deletions=0, graph_updates=0,
        physical_solver_evaluations=0, nnps_evaluations=0)
    if report['status'] != PASS:
        raise VisibleSupportNumericalDomainError(report)
    return report


def require_remote_filter_consistency(domain, graph):
    """Proved remote support must survive the unchanged actual endpoint filter."""
    if domain.get('classification') != DOMAIN_CLASSIFICATION or domain['status'] != PASS:
        raise ValueError('successful finite-tool endpoint domain report required')
    n = domain['edge_count']
    kept = np.asarray(graph['kept_raw_edge_indices'], dtype=int)
    blocked = np.asarray(graph['blocked_raw_edge_indices'], dtype=int)
    if (len(kept)+len(blocked) != n or set(kept) & set(blocked)
            or set(kept) | set(blocked) != set(range(n))):
        raise ValueError('actual filter does not partition the audited raw rows')
    if not set(domain['proven_visible_raw_edge_indices']) <= set(kept):
        raise ArithmeticError('unchanged filter discarded proved remote SPH support')


def static_path_adaptation_audit():
    """Only replace the query-domain call and assert its original query result."""
    expected_call = ast.parse('check_finite_partition_numerical_domain(start, direction, R, '
        'min(radius,a,height), half_bore=A, zlo=zlo, zhi=zhi)', mode='eval').body
    replacement = ast.parse('check_partition_numerical_domain(start,direction,R,'
        'min(radius,a,height))', mode='eval').body
    expected_assert = ast.parse('require_remote_partition_consistency(domain, result)').body[0]
    counts = dict(domain_calls=0, consistency_assertions=0)
    class Normalize(ast.NodeTransformer):
        def visit_Call(self, node):
            if ast.dump(node) == ast.dump(expected_call):
                counts['domain_calls'] += 1
                return replacement
            return self.generic_visit(node)
        def visit_Expr(self, node):
            if ast.dump(node) == ast.dump(expected_assert):
                counts['consistency_assertions'] += 1
                return None
            return self.generic_visit(node)
    actual = ast.parse(textwrap.dedent(inspect.getsource(audit_visible_support_path)))
    actual = Normalize().visit(actual)
    old = ast.parse(textwrap.dedent(inspect.getsource(legacy_path.audit_visible_support_path)))
    passed = (counts == dict(domain_calls=1, consistency_assertions=1)
        and ast.dump(actual) == ast.dump(old))
    if not passed:
        raise AssertionError('tight-path static adaptation exceeded registered domain-only changes')
    return dict(passed=True, **counts, unchanged_near_tool_guard=True,
        unchanged_path_bounds=True, unchanged_query_budget=True,
        unchanged_strict_geometry=True, edge_deletions=0,
        classification=CLASSIFICATION)


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
