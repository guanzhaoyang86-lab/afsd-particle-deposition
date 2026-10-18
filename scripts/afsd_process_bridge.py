"""Pure orchestration for a four-call real-geometry response measurement.

This supplies numerical trial perturbations, NEVER material trajectories.
No physical/backend functions are imported here.
"""
import copy
import hashlib
import numpy as np

FIELDS=('position_m','velocity_m_s','density_kg_m3','temperature_K','mean_velocity_m_s')
EPSILON=1e-7  # m/s; registered once, no epsilon search.
OLD_COST_SECONDS=4.04492670000036
LIMIT_SECONDS=50.
MATRIX=np.array([[1.,.25,.125],[.25,-.5,.0625],[.125,.0625,-.5]])
LABELS=('base','plus','minus','base_repeat')


def direction(position,mass,particle_id,free_axes):
    x=np.asarray(position,float);m=np.asarray(mass,float);ids=np.asarray(particle_id);free=np.asarray(free_axes,float)
    if m.ndim!=1:raise ValueError('mass vector required')
    n=len(m)
    if (n==0 or x.shape!=(n,3) or m.shape!=(n,) or ids.shape!=(n,) or free.shape!=(n,3)
        or not np.isfinite(x).all() or not np.isfinite(m).all() or not np.isfinite(ids).all()
        or not np.isfinite(free).all() or np.any(m<=0) or len(np.unique(ids))!=n
        or not np.all(ids==np.floor(ids)) or not np.all((free==0)|(free==1))):
        raise ValueError('finite unique particle IDs, positive masses, 3D positions and binary free axes required')
    order=np.argsort(ids,kind='stable')
    center=np.sum(m[order,None]*x[order],axis=0)/np.sum(m[order])
    length=float(np.max(np.max(x,axis=0)-np.min(x,axis=0)))
    if not np.isfinite(length) or length<=0:raise ValueError('nonzero physical extent required')
    raw=((x-center)/length)@MATRIX.T*free
    scale=float(np.max(np.linalg.norm(raw,axis=1)))
    if not np.isfinite(scale) or scale<=0:raise ValueError('nonzero free trial direction required')
    return raw/scale


def candidate_copy(candidate):
    if not isinstance(candidate,dict) or set(candidate)!=set(FIELDS):
        raise ValueError('exactly five candidate fields required')
    result={k:np.asarray(candidate[k],float).copy() for k in FIELDS}
    if result['density_kg_m3'].ndim!=1:raise ValueError('density vector required')
    n=len(result['density_kg_m3'])
    for key,value in result.items():
        shape=(n,) if key in ('density_kg_m3','temperature_K') else (n,3)
        if not n or value.shape!=shape or not np.isfinite(value).all():raise ValueError('invalid candidate '+key)
    return result


def candidate_hash(candidate):
    data=candidate_copy(candidate);digest=hashlib.sha256()
    for key in FIELDS:
        value=data[key];digest.update(key.encode());digest.update(str(value.shape).encode());digest.update(value.tobytes())
    return digest.hexdigest()


def perturb(candidate,D,dt,sign):
    c=candidate_copy(candidate);d=np.asarray(D,float)
    if (not np.isscalar(sign) or isinstance(sign,(bool,np.bool_)) or sign not in (-1,0,1)
        or not isinstance(dt,(float,int,np.floating,np.integer)) or isinstance(dt,(bool,np.bool_))
        or not np.isfinite(dt) or dt<=0
        or d.shape!=c['position_m'].shape or not np.isfinite(d).all()):
        raise ValueError('fixed finite direction, positive dt and sign -1/0/+1 required')
    if sign:
        increment=sign*EPSILON*d
        c['mean_velocity_m_s']+=increment
        c['velocity_m_s']+=2*increment
        c['position_m']+=dt*increment
    return c


def forecast(elapsed,remaining,measured):
    if (not np.isfinite(elapsed) or elapsed<0 or type(remaining) is not int or not 0<=remaining<=4
        or any(not np.isfinite(t) or t<0 for t in measured)):
        raise ValueError('valid elapsed, remaining count and measured times required')
    return float(elapsed+1.5*remaining*max([OLD_COST_SECONDS]+list(measured))+5.)


def may_continue(elapsed,remaining,measured,limit_seconds=LIMIT_SECONDS):
    if (not isinstance(limit_seconds,(float,int,np.floating,np.integer))
        or isinstance(limit_seconds,(bool,np.bool_)) or not np.isfinite(limit_seconds)
        or not 0<limit_seconds<=600):raise ValueError('explicit finite budget up to600s required')
    return forecast(elapsed,remaining,measured)<=limit_seconds


class BridgeFailure(RuntimeError):
    def __init__(self,index,label,candidate,records,original):
        super().__init__(str(original));self.index=index;self.label=label
        self.candidate=candidate_copy(candidate);self.records=copy.deepcopy(records);self.original=original


def run_four(evaluator,candidate,D,dt,clock,on_result,started=None,limit_seconds=LIMIT_SECONDS):
    """evaluator(C,dt), on_result(index,label,C,evaluation,seconds)->small receipt.

    Every C is freshly copied from the original, never the previous map. The
    caller persists full output and verifies accepted-state purity in on_result.
    A callback/evaluation failure is terminal; there is no retry or extra call.
    """
    base=candidate_copy(candidate);base_hash=candidate_hash(base)
    first=clock() if started is None else started
    measured=[];records=[]
    for index,(label,sign) in enumerate(zip(LABELS,(0,1,-1,0)),1):
        elapsed=clock()-first
        if not may_continue(elapsed,5-index,measured,limit_seconds):
            return dict(status='PARTIAL_BUDGET_NOT_RUN_REMAINDER',records=records,
                completed=len(records),not_run=list(LABELS[index-1:]),
                next_forecast_seconds=forecast(elapsed,5-index,measured))
        c=perturb(base,D,dt,sign);before=candidate_hash(c);tick=clock()
        try:
            evaluation=evaluator(c,dt);seconds=clock()-tick
            if seconds<0 or not np.isfinite(seconds):raise ValueError('invalid monotonic measurement')
            if candidate_hash(c)!=before:raise RuntimeError('evaluator mutated trial input')
            receipt=on_result(index,label,c,evaluation,seconds)
            if candidate_hash(c)!=before or candidate_hash(base)!=base_hash:
                raise RuntimeError('callback mutated trial input')
            full_seconds=clock()-tick
            if not np.isfinite(full_seconds) or full_seconds<seconds:
                raise ValueError('invalid complete evaluation/archive measurement')
            measured.append(full_seconds)
            records.append(dict(index=index,label=label,evaluation_seconds=seconds,
                complete_measurement_seconds=full_seconds,receipt=receipt))
        except Exception as error:
            raise BridgeFailure(index,label,c,records,error) from error
    return dict(status='FOUR_EVALUATIONS_COMPLETED_NOT_ACCEPTANCE',records=records,
        completed=4,not_run=[],elapsed_seconds=clock()-first)
