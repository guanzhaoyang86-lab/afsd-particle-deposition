"""Independent contact-integral limits; no backend, particles or trajectory.

The folded force is a declared mathematical PRIVATE extension. Its integral
is not assumed differentiable merely because every fixed Gauss node is.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import threading
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT/'docs/PROCESS_RESPONSE_LIMITS_CONTRACT.md'
FREEZE = ROOT/'docs/PROCESS_STOP_FREEZE.json'
ORDERS = (8, 16, 32, 64, 128)
_RULES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, np.ndarray):return value.tolist()
    if isinstance(value, np.generic):return value.item()
    if isinstance(value, dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):return [clean(v) for v in value]
    return value


def write(path, value):
    path=Path(path);pending=path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    pending.replace(path)


def verify_frozen():
    manifest=json.loads(FREEZE.read_text(encoding='utf-8'))
    records=manifest['files_sha256']
    if len(records)!=4516:raise ValueError('registered4516-file frozen inventory differs')
    changed=[]
    for name,digest in records.items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file() or sha(path)!=digest:
            changed.append(name)
    if changed:raise ValueError('frozen source/evidence changed: '+repr(changed))
    return dict(files_checked=len(records),changed_files=changed,manifest_sha256=sha(FREEZE))


def rule(order):
    if order not in _RULES:
        x,w=np.polynomial.legendre.leggauss(order)
        _RULES[order]=(.5*(x+1.),.5*w)
    return _RULES[order]


def closed_impulse(epsilon, s0=.4, C=1., h=1., U=1., folded=True):
    """Direct elementary antiderivatives, with explicit measure-zero limit."""
    a=float(s0);b=1.-a;rho=float(epsilon)/U
    if not 0<=a<=1 or min(C,h,U)<=0:raise ValueError('registered positive SI inputs')
    if rho==0.:return np.array([-C*h if folded else -C*h*(b-a),0.])
    radius=abs(rho)
    if folded:
        return -C*h*np.array([math.hypot(a,rho)+math.hypot(b,rho)-2*radius,
            rho*(math.asinh(b/radius)-math.asinh(a/radius))])
    return -C*h*np.array([math.hypot(b,rho)-math.hypot(a,rho),
        rho*(math.asinh(b/radius)+math.asinh(a/radius))])


def stable_parallel_delta(epsilon,s0=.4,C=1.,h=1.,U=1.):
    """Jx(epsilon)-Jx(0), avoiding cancellation; not heat clipping."""
    rho=float(epsilon)/U
    if rho==0.:return 0.
    a=float(s0);b=1.-a
    return C*h*(2*abs(rho)-rho*rho/(math.hypot(a,rho)+a)
                -rho*rho/(math.hypot(b,rho)+b))


def transformed_direct_force_quadrature(epsilon,s0,order,C=1.,h=1.,U=1.):
    """Independent direct force samples after q=+/-|rho|sinh(z).

    The closed impulse routine is not called here. This substitution resolves
    the small closest-approach region without adding a slip floor or clipping.
    """
    if epsilon==0.:raise ValueError('nonzero epsilon required by transformed rule')
    rho=epsilon/U;radius=abs(rho);nodes,weights=rule(order);answer=np.zeros(2)
    for sign,extent in ((-1.,s0),(1.,1.-s0)):
        if extent==0.:continue
        end=math.asinh(extent/radius);z=end*nodes
        q=sign*radius*np.sinh(z)
        velocity=np.column_stack((U*q,np.full(len(q),epsilon)))
        force=-C*np.sign(q)[:,None]*velocity/np.linalg.norm(velocity,axis=1)[:,None]
        ds_dz=radius*np.cosh(z)
        answer+=h*end*np.sum((weights*ds_dz)[:,None]*force,axis=0)
    return answer


def fixed_gauss_impulse(epsilon,s0,order,C=1.,h=1.,U=1.):
    nodes,weights=rule(order);q=nodes-s0
    if np.any(q==0.):raise ValueError('a Gauss node is on the singular branch point')
    velocity=np.column_stack((U*q,np.full(len(q),epsilon)))
    force=-C*np.sign(q)[:,None]*velocity/np.linalg.norm(velocity,axis=1)[:,None]
    return h*np.sum(weights[:,None]*force,axis=0)


def fixed_gauss_parallel_delta(epsilon,s0,order,C=1.,h=1.,U=1.):
    nodes,weights=rule(order);q=np.abs(nodes-s0);rho=epsilon/U
    radius=np.hypot(q,rho)
    return C*h*float(np.sum(weights*rho*rho/(radius*(radius+q))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh workspace output required')
    out.mkdir(parents=True);started=time.perf_counter()
    sources={Path(__file__).relative_to(ROOT).as_posix():sha(__file__),
        CONTRACT.relative_to(ROOT).as_posix():sha(CONTRACT)}
    report=dict(status='RUNNING_MATHEMATICAL_LIMIT_PROBE',source_sha256=sources,
        classification='INDEPENDENT_PRESCRIBED_PATH_MATH_NOT_MATERIAL_TRAJECTORY',
        new_afsd_steps=0,new_particle_steps=0,accepted=False,restart_allowed=False,
        production_authorized=False,checks=[],quadrature_checks=[],interior_limits=[],
        endpoint_limits=[],finite_gauss_maps=[],ordinary_law_limits=[])
    def timeout():
        write(out/'HARD_TIMEOUT.json',dict(status='TIMEOUT55S',accepted=False));os._exit(124)
    timer=threading.Timer(55.,timeout);timer.daemon=True;timer.start();code=2
    def check(name,error,bound):
        passed=bool(np.isfinite(error) and error<=bound)
        report['checks'].append(dict(name=name,error=float(error),bound=float(bound),passed=passed))
        if not passed:raise ValueError('predeclared mathematical check failed: '+name)
    try:
        report['frozen_before']=verify_frozen()
        for s0 in (.4,.5,0.,1.):
            for magnitude in (.2,.02,.002,.0002,1e-8):
                for epsilon in (magnitude,-magnitude):
                    exact=closed_impulse(epsilon,s0)
                    q96=transformed_direct_force_quadrature(epsilon,s0,96)
                    q192=transformed_direct_force_quadrature(epsilon,s0,192)
                    check('closed_vs_direct_%s_%s'%(s0,epsilon),float(np.max(np.abs(q192-exact))),2e-12)
                    check('direct96_vs192_%s_%s'%(s0,epsilon),float(np.max(np.abs(q192-q96))),2e-12)
                    report['quadrature_checks'].append(dict(s0=s0,epsilon_m_s=epsilon,
                        closed_impulse_Ns=exact,direct96_Ns=q96,direct192_Ns=q192))
        for magnitude in (1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14):
            positive=stable_parallel_delta(magnitude)/magnitude
            negative=stable_parallel_delta(-magnitude)/(-magnitude)
            transverse=closed_impulse(magnitude)[1]/magnitude
            report['interior_limits'].append(dict(epsilon_m_s=magnitude,
                right_parallel_quotient=positive,left_parallel_quotient=negative,
                centered_parallel_quotient=0.,transverse_quotient=transverse,
                reference_right=2.,reference_left=-2.,reference_transverse=-math.log(1.5)))
            if magnitude<=1e-8:
                check('interior_right_'+str(magnitude),abs(positive-2.),3e-8)
                check('interior_left_'+str(magnitude),abs(negative+2.),3e-8)
                check('interior_transverse_'+str(magnitude),abs(transverse+math.log(1.5)),1e-11)
            for s0 in (0.,1.):
                impulse=closed_impulse(magnitude,s0)
                slope=stable_parallel_delta(magnitude,s0)/magnitude
                sign=-1. if s0==0. else 1.
                asymptote=sign*math.log(2./magnitude)
                report['endpoint_limits'].append(dict(s0=s0,epsilon_m_s=magnitude,
                    right_parallel_quotient=slope,left_parallel_quotient=-slope,
                    transverse_quotient=impulse[1]/magnitude,logarithmic_asymptote=asymptote))
                if magnitude<=1e-8:
                    check('endpoint_parallel_%s_%s'%(s0,magnitude),abs(slope-1.),1e-8)
                    check('endpoint_log_%s_%s'%(s0,magnitude),abs(impulse[1]/magnitude-asymptote),1e-11)
            ordinary=closed_impulse(magnitude,.4,folded=False)
            report['ordinary_law_limits'].append(dict(epsilon_m_s=magnitude,
                transverse_quotient=ordinary[1]/magnitude,
                leading_logarithmic_asymptote=-math.log(4*.4*.6/(magnitude*magnitude))))
        for s0 in (.4,.5):
            for order in ORDERS:
                nodes,weights=rule(order);q=nodes-s0
                epsilon=1e-4*float(np.min(np.abs(q)))
                exact_discrete=np.array([0.,-float(np.sum(weights/q))])
                fd=(fixed_gauss_impulse(epsilon,s0,order)-fixed_gauss_impulse(-epsilon,s0,order))/(2*epsilon)
                for axis in range(2):
                    check('fixed_rule_fd_%s_%s_%s'%(s0,order,axis),abs(fd[axis]-exact_discrete[axis]),
                        3e-7*max(1.,abs(exact_discrete[axis])))
                report['finite_gauss_maps'].append(dict(s0=s0,order=order,fd_epsilon_m_s=epsilon,
                    exact_finite_rule_derivative=exact_discrete,centered_finite_difference=fd,
                    continuous_transverse_derivative=-math.log((1-s0)/s0),
                    continuous_full_vector_derivative='NOT_DIFFERENTIABLE'))
        finite=fixed_gauss_parallel_delta(1e-10,.4,32)/1e-10
        continuous=stable_parallel_delta(1e-10,.4)/1e-10
        check('fixed32_parallel_limit_below1e-6',finite,1e-6)
        check('continuous_parallel_limit_above1.99',max(0.,1.99-continuous),0.)
        report['noncommuting_limits']=dict(order=32,epsilon_m_s=1e-10,
            finite_rule_parallel_quotient=finite,continuous_parallel_quotient=continuous,
            finite_rule_classical_parallel_derivative=0.,continuous_right_parallel_derivative=2.,
            continuous_left_parallel_derivative=-2.,centered_difference_is_not_cusp_test=True)
        report['frozen_after']=verify_frozen()
        if any(sha(ROOT/name)!=digest for name,digest in sources.items()):raise ValueError('probe source drift')
        report.update(status='PASS_MATHEMATICAL_NONDIFFERENTIABILITY_EVIDENCE_NOT_SOLVER_VALIDATION',
            continuous_cross_zero_vector_jacobian='NOT_DIFFERENTIABLE',
            classical_continuous_cross_zero_newton_claim='BLOCKED',
            named_fixed_rule_preconditioner='CONDITIONALLY_ELIGIBLE_FOR_SEPARATELY_REGISTERED_PRIVATE_TESTS',
            saved_trial005_diagnostic='NOT_RUN',trajectory_status='NOT_RUN');code=0
    except Exception as error:
        report.update(status='FAIL_LIMIT_PROBE_OR_INTEGRITY',error=repr(error))
    finally:
        report['wall_seconds']=time.perf_counter()-started
        report['checks_run']=len(report['checks'])
        write(out/'report.json',report)
        timer.cancel()
    print(json.dumps(clean({k:v for k,v in report.items() if k not in (
        'checks','quadrature_checks','interior_limits','endpoint_limits','finite_gauss_maps',
        'ordinary_law_limits','source_sha256')})),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())
