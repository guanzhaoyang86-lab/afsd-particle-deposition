"""Explicit engineering numerical domain for the inherited solid partitioner.

The frozen partitioner forms bb**2 - 4*aa*cc. Stable root extraction does not
undo cancellation already present in this discriminant. This independent
guard rejects unresolved signs, and limits coordinates, segment length, and
the queried cylinder radius to eight times the original smallest tool scale.
It does not modify edges, dimensions, or the old frozen filter.

The error envelope below includes conservative coefficient and discriminant
roundoff allowances. This is a declared floating-point engineering domain and
adversarially tested guard, not a formal directed-rounding interval proof.
"""
import numpy as np
from visible_support import _inputs

PASS='PASS_VISIBLE_SUPPORT_NUMERICAL_DOMAIN'
UNRESOLVED='UNRESOLVED_VISIBLE_SUPPORT_NUMERICAL_DOMAIN'
MAXIMUM_SCALE_RATIO=8.


class VisibleSupportNumericalDomainError(ArithmeticError):
    def __init__(self,report):
        self.report=report
        super().__init__('visible support numerical domain unresolved: '+','.join(report.get('reasons',[])))


def check_partition_numerical_domain(start,direction,radius,minimum_tool_dimension):
    """Check the exact arrays/radius about to reach a partition query.

    ``minimum_tool_dimension`` is the ORIGINAL physical tool scale, even when
    ``radius`` belongs to an expanded enclosure. No geometric query is made.
    """
    s=np.asarray(start,float);d=np.asarray(direction,float)
    radius=float(radius);scale=float(minimum_tool_dimension)
    if s.ndim!=2 or s.shape[1]!=3 or d.shape!=s.shape or not np.isfinite(s).all() or not np.isfinite(d).all():
        raise ValueError('finite matching N by 3 segment arrays required')
    if not np.isfinite([radius,scale]).all() or min(radius,scale)<=0:
        raise ValueError('positive finite queried radius and original tool scale required')
    report=dict(status=PASS,maximum_allowed_scale_ratio=MAXIMUM_SCALE_RATIO,
        original_minimum_tool_dimension_m=scale,queried_cylinder_radius_m=radius,edge_count=len(s),
        reasons=[],uncertain_edge_indices=[],physical_solver_evaluations=0,
        limitation='Engineering floating-point error envelope; not formal interval arithmetic.')
    with np.errstate(over='ignore',invalid='ignore',under='ignore',divide='ignore'):
        ratio=max(float(np.linalg.norm(s,axis=1).max(initial=0.))/scale,
            float(np.linalg.norm(s+d,axis=1).max(initial=0.))/scale,
            float(np.linalg.norm(d,axis=1).max(initial=0.))/scale,radius/scale)
    report['maximum_actual_scale_ratio']=ratio if np.isfinite(ratio) else None
    if not np.isfinite(ratio) or ratio>MAXIMUM_SCALE_RATIO:
        report['status']=UNRESOLVED;report['reasons'].append('COORDINATE_EDGE_OR_RADIUS_SCALE_RATIO_EXCEEDS_EIGHT')
        return report
    eps=np.finfo(float).eps
    with np.errstate(over='ignore',invalid='ignore',under='ignore',divide='ignore'):
        aa=np.sum(d[:,:2]**2,axis=1)
        bb=2*np.sum(s[:,:2]*d[:,:2],axis=1)
        squared_start=np.sum(s[:,:2]**2,axis=1)
        cc=squared_start-radius**2
        disc=bb*bb-4*aa*cc
        ea=8*eps*aa
        eb=16*eps*np.sum(np.abs(s[:,:2]*d[:,:2]),axis=1)
        ec=8*eps*(squared_start+radius**2)
        bound=8*eps*(bb*bb+np.abs(4*aa*cc))+(2*np.abs(bb)+eb)*eb+4*((np.abs(aa)+ea)*ec+np.abs(cc)*ea)
        bound=np.nextafter(bound,np.inf)
    finite=np.isfinite(aa)&np.isfinite(bb)&np.isfinite(cc)&np.isfinite(disc)&np.isfinite(bound)
    # aa=0 denotes an axis-parallel edge and no radial roots are needed. Its
    # constant radial classification must still be separated from the boundary.
    uncertain=(~finite)|((aa>0)&(np.abs(disc)<=bound))|((aa==0)&(np.abs(cc)<=ec))
    indices=np.flatnonzero(uncertain)
    report['uncertain_edge_indices']=indices.astype(int).tolist()
    report['cylinder_discriminant_uncertain_count']=len(indices)
    if len(indices):
        report['status']=UNRESOLVED;report['reasons'].append('CYLINDER_SIGN_UNRESOLVED_BY_ROUNDOFF_ENVELOPE')
    return report


def require_visible_support_numerical_domain(x,i,j,pose,p,particle_ids):
    """Run before the frozen endpoint filter; raise on ANY unresolved domain.

    Both strict and actual roundoff-expanded partition inputs are checked.
    The caller must reject the entire candidate on failure, not drop extra rows.
    """
    _,i,j,_,_,_,local,radius,a,height,margin=_inputs(x,i,j,pose,p,particle_ids)
    start=local[i];direction=local[j]-start;scale=min(radius,a,height)
    strict=check_partition_numerical_domain(start,direction,radius,scale)
    expanded=check_partition_numerical_domain(start,direction,np.nextafter(radius+margin,np.inf),scale)
    report=dict(status=PASS if strict['status']==PASS and expanded['status']==PASS else UNRESOLVED,
        scope='STRICT_AND_ROUNDOFF_EXPANDED_ENDPOINT_PARTITION_INPUTS',
        strict=strict,expanded=expanded,reasons=sorted(set(strict['reasons']+expanded['reasons'])))
    if report['status']!=PASS:raise VisibleSupportNumericalDomainError(report)
    return report
