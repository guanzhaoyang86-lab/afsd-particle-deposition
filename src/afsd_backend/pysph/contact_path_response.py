"""Fixed-quadrature numerical contact response, never a force/heat model.

The smooth point formulas are valid away from full zero and branch switches.
The finite-rule derivative is NOT a certificate of differentiability of the
continuous cross-zero impulse. In particular, PRIVATE branch continuation
has a transverse cusp. This component cannot authorize an event or trajectory.
"""
import numpy as np
from .hht_path import parameters


def _scalar(value,name,positive=False):
    if isinstance(value,(bool,np.bool_)) or np.ndim(value)!=0:
        raise ValueError('finite scalar '+name+' required')
    value=float(value)
    if not np.isfinite(value) or (value<=0 if positive else value<0):
        raise ValueError('invalid '+name)
    return value


def _vector(value,name):
    result=np.array(value,float,copy=True)
    if result.shape!=(3,) or not np.isfinite(result).all():raise ValueError('finite Cartesian '+name)
    return result


def point_response(normal,slip_m_s,normal_force_N,normal_stiffness_N_m,friction,branch_sign=1):
    """Derivatives dF/dx and dF/dv of one fixed-normal contact, SI units.

    Fn=kn*gap, gap derivative wrt material x is -n. sigma=+1 means physical
    Ft=-mu*Fn*u/|u|; sigma=-1 is ONLY the private extended-side equation.
    Neither a force nor an updated physical field is returned.
    """
    n=_vector(normal,'normal');u=_vector(slip_m_s,'tangential slip')
    fn=_scalar(normal_force_N,'normal force',True);kn=_scalar(normal_stiffness_N_m,'normal stiffness',True)
    mu=_scalar(friction,'friction')
    if isinstance(branch_sign,(bool,np.bool_)) or np.ndim(branch_sign)!=0 or branch_sign not in (-1,1):
        raise ValueError('branch_sign must be +1 physical or -1 private extension')
    eps=np.finfo(float).eps;speed=float(np.linalg.norm(u))
    if abs(np.linalg.norm(n)-1)>512*eps:raise ValueError('unit normal required')
    if not np.isfinite(speed) or speed==0:raise ValueError('resolved nonzero slip required; no zero regularization')
    if abs(n@u)>512*eps*max(speed,np.finfo(float).tiny):raise ValueError('full tangential slip required')
    tangent=u/speed;projector=np.eye(3)-np.outer(n,n)
    cap=mu*fn
    kx=-kn*np.outer(n,n)+branch_sign*mu*kn*np.outer(tangent,n)
    kv=-branch_sign*(cap/speed)*(projector-np.outer(tangent,tangent))
    if not np.isfinite(cap) or not np.isfinite(kx).all() or not np.isfinite(kv).all():
        raise ValueError('nonfinite contact response; no silent cap')
    return dict(Kx_N_m=kx,Kv_N_s_m=kv,slip_speed_m_s=speed,branch_sign=int(branch_sign),
        classification='SMOOTH_POINT_CONTACT_RESPONSE_NOT_PHYSICAL_STATE',
        accepted=False,restart_allowed=False,continuous_path_differentiability_certified=False)


def path_response(nodes,end_row,dt,weight=.9):
    """Full DeltaJ contribution to a FORCE-unit HHT root, fixed-rule derivative.

    Rcontact=-a*F1-(1-a)*F0-(Jactual-h*(F0+F1)/2)/(gamma*h).
    Left state and prescribed wall programme are fixed. The material/NNPS/
    density/thermal derivative is NOT included. Each result is one3x3 block.
    Nodes have normalized weights on[0,1], not time-weighted impulses.
    """
    h=_scalar(dt,'dt',True);gamma,beta=parameters(weight);a=float(weight);d=beta-gamma/2
    rows=list(nodes)
    if not rows:raise ValueError('nonempty normalized quadrature required')
    names=('normal','tangential_relative_velocity_m_s','normal_force_N',
           'normal_stiffness_N_m','friction','branch_sign')
    def point(row):
        if any(key not in row for key in names):raise ValueError('complete source-bound quadrature row required')
        return point_response(row['normal'],row['tangential_relative_velocity_m_s'],
            row['normal_force_N'],row['normal_stiffness_N_m'],row['friction'],row['branch_sign'])
    end=point(end_row);total=np.zeros((3,3));weights=[];signs=[];node_blocks=[]
    for row in rows:
        if not np.array_equal(np.asarray(row['normal']),np.asarray(end_row['normal'])):
            raise ValueError('fixed normal required throughout the interval')
        if row['normal_stiffness_N_m']!=end_row['normal_stiffness_N_m'] or row['friction']!=end_row['friction']:
            raise ValueError('fixed stiffness/friction required')
        s=_scalar(row['fraction'],'fraction');w=_scalar(row['weight'],'quadrature weight',True)
        if s>1:raise ValueError('quadrature fraction outside[0,1]')
        response=point(row)
        va=h*(gamma*s+6*d*s*(1-s))
        xa=h*h*((gamma/2+3*d)*s*s-2*d*s*s*s)
        block=response['Kx_N_m']*xa+response['Kv_N_s_m']*va
        total+=h*w*block;weights.append(w);signs.append(row['branch_sign']);node_blocks.append(block)
    if abs(sum(weights)-1)>1e-12:raise ValueError('quadrature weights must sum to1')
    endpoint=end['Kx_N_m']*(beta*h*h)+end['Kv_N_s_m']*(gamma*h)
    root=(1/(2*gamma)-a)*endpoint-total/(gamma*h)
    if not np.isfinite(root).all():raise ValueError('nonfinite full contact residual response')
    return dict(actual_impulse_acceleration_response_kg_s=total,
        endpoint_force_acceleration_response_kg=endpoint,
        root_contact_response_kg=root,node_force_acceleration_response_kg=np.asarray(node_blocks),
        classification='FROZEN_QUADRATURE_NUMERICAL_RESPONSE_NOT_CONTINUOUS_CROSS_ZERO_TANGENT',
        impulse_response_units='kg*s',force_response_units='kg',
        accepted=False,restart_allowed=False,physical_force_changed=False,
        continuous_path_differentiability_certified=False,
        quadrature_selection_derivative_included=False,
        multiple_sampled_branch_signs=len(set(signs+[end_row['branch_sign']]))>1)
