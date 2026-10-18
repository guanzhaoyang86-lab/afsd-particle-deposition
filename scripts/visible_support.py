"""Pure endpoint visibility filter for one finite open-square cylindrical tool.

This is a NEW spatial graph policy, not the old 1 nm diagnostic guard. Every
strict positive-length solid intersection is blocked. A second intersection
with a roundoff-expanded solid also blocks uncertain grazing edges. Boundary
uncertainty therefore removes support; it never authorizes a through-solid
edge. The existing checked analytic partition implementation is reused.

Input edges must already have canonical undirected orientation i<j. Their raw
order is preserved, so a completely visible graph and its supplied gradients
are returned bit-for-bit unchanged. No NNPS, mechanics or heat is evaluated.
"""
import hashlib
import json
import numpy as np
from audit_solid_occluded_neighbors import _partitions

CLASSIFICATION='ENGINEERING_APPROXIMATE_ENDPOINT_VISIBLE_SUPPORT_V1'
POLICY='BLOCK_ANY_POSITIVE_STRICT_INTERSECTION_OR_ROUNDOFF_EXPANDED_SOLID_HIT'


def _hash_array(value):
    a=np.ascontiguousarray(value)
    h=hashlib.sha256()
    h.update(str(a.dtype).encode('ascii'));h.update(str(a.shape).encode('ascii'));h.update(a.tobytes())
    return h.hexdigest()


def _topology_hash(pairs):
    pairs=np.sort(np.asarray(pairs,dtype='<i8').reshape(-1,2),axis=1)
    if len(pairs):pairs=pairs[np.lexsort((pairs[:,1],pairs[:,0]))]
    return _hash_array(pairs)


def _inputs(x,i,j,pose,p,particle_ids):
    x=np.asarray(x,dtype=float);i=np.asarray(i);j=np.asarray(j);ids=np.asarray(particle_ids)
    if x.ndim!=2 or x.shape[1]!=3 or not len(x) or not np.isfinite(x).all():
        raise ValueError('finite nonempty N by 3 positions required')
    if i.ndim!=1 or j.shape!=i.shape or i.dtype.kind not in 'iu' or j.dtype.kind not in 'iu':
        raise ValueError('matching one-dimensional integer edge indices required')
    if np.any(i<0) or np.any(j>=len(x)) or np.any(i>=j):
        raise ValueError('canonical edges require 0 <= i < j < particle count')
    if len(i) and len(np.unique(np.column_stack((i,j)),axis=0))!=len(i):
        raise ValueError('duplicate canonical support edge')
    if (ids.shape!=(len(x),) or ids.dtype.kind not in 'iu' or np.any(ids<0)
            or np.any(ids>np.iinfo(np.int64).max) or len(np.unique(ids))!=len(ids)):
        raise ValueError('unique nonnegative persistent integer IDs required')
    origin=np.asarray(pose['origin'] if isinstance(pose,dict) else pose.origin,dtype=float)
    rotation=np.asarray(pose['rotation'] if isinstance(pose,dict) else pose.rotation,dtype=float)
    if origin.shape!=(3,) or rotation.shape!=(3,3) or not np.isfinite(origin).all() or not np.isfinite(rotation).all():
        raise ValueError('finite origin and rotation required')
    defect=max(float(np.max(np.abs(rotation.T@rotation-np.eye(3)))),abs(float(np.linalg.det(rotation))-1.))
    if defect>1e-12:raise ValueError('proper orthogonal tool rotation required')
    radius=float(p['tool_radius']);half_bore=.5*float(p['bore_width']);height=float(p['tool_top'])-float(p['layer'])
    if not np.isfinite([radius,half_bore,height]).all() or min(radius,half_bore,height)<=0:
        raise ValueError('positive finite tool dimensions required')
    relative=x-origin
    local=relative@rotation
    scale=max(float(np.max(np.abs(x))),float(np.max(np.abs(origin))),radius,half_bore,height,
              float(np.max(np.abs(relative))),float(np.max(np.abs(local))))
    # Covers subtraction, 3-term dot products and analytic intersection tests.
    # Absolute world-coordinate scale deliberately exposes cancellation from a
    # large translation; orthogonality error is separately enlarged.
    margin=256*np.finfo(float).eps*scale+8*defect*float(np.max(np.linalg.norm(relative,axis=1)))
    margin=float(np.nextafter(margin,np.inf))
    if not np.isfinite(local).all() or not np.isfinite(margin) or margin>=min(radius,half_bore,height)/16:
        raise ArithmeticError('tool geometry unresolved at supplied coordinate/rotation scale')
    if np.any(np.linalg.norm(local[j]-local[i],axis=1)==0):
        raise ValueError('coincident support endpoints are unsupported')
    return x,i,j,ids,origin,rotation,local,radius,half_bore,height,margin


def filter_visible_support(x,i,j,g,pose,p,particle_ids):
    """Return a conservative i/j/g subset and immutable-identity diagnostics.

    Tags and component labels are deliberately absent from this API. Persistent
    IDs affect hashes/diagnostics only and do not influence the geometric mask.
    Empty edge lists are allowed. Inputs are never modified. The output proves
    only endpoint geometry, never visibility over the preceding material path.
    """
    with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
        x,i,j,ids,origin,rotation,local,radius,half_bore,height,margin=_inputs(x,i,j,pose,p,particle_ids)
        gradient=np.asarray(g)
        if gradient.shape!=(len(i),3) or gradient.dtype.kind!='f' or not np.isfinite(gradient).all():
            raise ValueError('finite floating edge gradients with shape E by 3 required')
        start=local[i];direction=local[j]-start
        cuts,inside,maximum,total=_partitions(start,direction,radius,half_bore,0.,height)
        _,_,expanded_maximum,expanded_total=_partitions(start,direction,
            np.nextafter(radius+margin,np.inf),np.nextafter(half_bore-margin,-np.inf),
            np.nextafter(-margin,-np.inf),np.nextafter(height+margin,np.inf))
    if not all(np.isfinite(a).all() for a in (maximum,total,expanded_maximum,expanded_total)):
        raise ArithmeticError('nonfinite analytic visibility result')
    strict=total>0.;expanded=expanded_total>0.;blocked=strict|expanded
    kept_index=np.flatnonzero(~blocked);blocked_index=np.flatnonzero(blocked)
    pid_pairs=np.column_stack((ids[i],ids[j])).astype('<i8',copy=False)
    records=[]
    for k in blocked_index:
        intervals=[[float(cuts[k,n]),float(cuts[k,n+1])] for n in range(inside.shape[1])
                   if inside[k,n] and cuts[k,n+1]>cuts[k,n]]
        records.append(dict(raw_edge_index=int(k),i=int(i[k]),j=int(j[k]),
            particle_id_i=int(ids[i[k]]),particle_id_j=int(ids[j[k]]),
            reason='STRICT_POSITIVE_LENGTH_SOLID_INTERSECTION' if strict[k] else 'ROUNDOFF_GRAZING_UNCERTAIN_CONSERVATIVELY_BLOCKED',
            strict_maximum_contiguous_length_m=float(maximum[k]),strict_total_length_m=float(total[k]),
            roundoff_expanded_maximum_contiguous_length_m=float(expanded_maximum[k]),
            strict_intervals=intervals))
    fingerprints=dict(particle_id_order_sha256=_hash_array(ids.astype('<i8')),
        raw_edge_pid_order_sha256=_hash_array(pid_pairs),raw_topology_sha256=_topology_hash(pid_pairs),
        visible_topology_sha256=_topology_hash(pid_pairs[kept_index]),blocked_topology_sha256=_topology_hash(pid_pairs[blocked_index]),
        raw_gradient_sha256=_hash_array(gradient),visible_gradient_sha256=_hash_array(gradient[kept_index]),
        world_positions_sha256=_hash_array(x),tool_origin_sha256=_hash_array(origin),tool_rotation_sha256=_hash_array(rotation))
    geometry_json=json.dumps(dict(radius=radius,half_bore=half_bore,height=height,policy=POLICY),sort_keys=True,separators=(',',':'))
    fingerprints['geometry_policy_sha256']=hashlib.sha256(geometry_json.encode('utf-8')).hexdigest()
    diagnostics=dict(classification=CLASSIFICATION,policy=POLICY,scope='ENDPOINT_CANONICAL_SUPPORT_GRAPH_ONLY',
        raw_edge_count=len(i),visible_edge_count=len(kept_index),blocked_edge_count=len(blocked_index),
        strict_blocked_edge_count=int(strict.sum()),roundoff_only_blocked_edge_count=int((blocked&~strict).sum()),
        physical_solid_intersection_length_tolerance_m=0.,old_1nm_guard_not_used_for_filter=True,
        roundoff_solid_expansion_m=margin,maximum_strict_intersection_length_m=float(maximum.max(initial=0.)),
        canonical_orientation_preserved=True,raw_edge_order_preserved=True,nnps_evaluations=0,physical_solver_evaluations=0,
        blocked_edges=records,limitations=['Endpoint filter only; a separate path certificate is needed for an interval.',
            'Exact grazing may be conservatively blocked by the explicit arithmetic expansion.',
            'No constitutive, mechanical, thermal, or continuum accuracy is established by visibility.'])
    return dict(i=i[kept_index].copy(),j=j[kept_index].copy(),g=gradient[kept_index].copy(),
        kept_raw_edge_indices=kept_index,blocked_raw_edge_indices=blocked_index,
        strict_blocked_raw_edge_indices=np.flatnonzero(strict),roundoff_blocked_raw_edge_indices=np.flatnonzero(blocked&~strict),
        raw_edge_pid_pairs=pid_pairs.copy(),visible_edge_pid_pairs=pid_pairs[kept_index].copy(),
        blocked_edge_pid_pairs=pid_pairs[blocked_index].copy(),diagnostics=diagnostics,fingerprints=fingerprints)
