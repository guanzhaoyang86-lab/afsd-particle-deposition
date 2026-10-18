"""Fixed-normal Q7 onset eligibility; unchanged J2/cone roots/atomic commit."""
import json
from pathlib import Path
import sys
from unittest.mock import patch
import numpy as np

import probe_process_plastic_friction as parent
from probe_process_plastic_friction_recovery import RecoveryBlock
from probe_process_friction_cone import feasible_cone
from probe_process_plastic import CandidateFailure

ROOT = parent.ROOT
CONFIG = ROOT/'configs/benchmarks/process_plastic_friction_onset.json'
EXTRA_SOURCES = ['scripts/probe_process_plastic_friction_onset.py',
    'configs/benchmarks/process_plastic_friction_onset.json',
    'scripts/probe_process_plastic_friction_recovery.py',
    'configs/benchmarks/process_plastic_friction_q7_recovery.json']


def boundary_eligibility(owner):
    N = np.asarray(owner.normal_force)
    tau = np.asarray(owner.tangent_force)
    if N.shape != (len(owner.nodes),) or tau.shape != (len(N),2):
        raise ValueError('invalid accepted contact shape')
    if not np.isfinite(N).all() or not np.isfinite(tau).all() or np.any(N<0):
        raise ValueError('invalid accepted contact state')
    capacity = owner.mu*N
    margin = capacity-np.linalg.norm(tau,axis=1)
    ft = owner.limits['force_tolerance_N']
    if np.any(margin < -ft):
        raise CandidateFailure('accepted state outside cone tolerance')
    return (capacity>ft) & (np.abs(margin)<=ft), margin


def advance_onset(owner, t, on_accept=None, trace=None):
    if not np.isfinite(t) or t<=owner.time:
        raise ValueError('time must advance')
    emit = trace if trace is not None else lambda item: None
    lim = owner.limits

    def branch(time, mask, stage):
        try:
            c = owner.branch(time,mask)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            emit(dict(stage=stage,time_s=time,released=mask.tolist(),status='BLOCKED',
                      error=str(exc),iterations=getattr(exc,'iterations',[])))
            raise
        emit(dict(stage=stage,time_s=time,released=mask.tolist(),status='ROOT_RETURNED_NOT_ACCEPTED',
            cone_margin_N=(owner.mu*c['normal_force']-np.linalg.norm(c['tangent_force'],axis=1)).tolist(),
            relative_increment=c['relative_increment'].tolist(),iterations=c.get('iteration_trace',[])))
        return c

    for _ in range(lim['maximum_events_per_interval']):
        eligible, oldmargin = boundary_eligibility(owner)
        probe = owner.time+(t-owner.time)*lim['mode_probe_fraction']
        if not owner.time<probe<=t:
            raise CandidateFailure('unresolved onset probe time')
        start = branch(probe,eligible,'BOUNDARY_PROBE')
        margin = owner.mu*start['normal_force']-np.linalg.norm(start['tangent_force'],axis=1)
        inward = eligible & (margin>lim['force_tolerance_N'])
        released = eligible & ~inward
        emit(dict(stage='CLASSIFICATION',base_time_s=owner.time,eligible=eligible.tolist(),
            initial_margin_N=oldmargin.tolist(),resolved_inward=inward.tolist(),released=released.tolist()))
        c = branch(t,released,'TARGET')
        if feasible_cone(c,released,owner.mu,lim):
            c.update(event_kind='REQUESTED')
        else:
            lo,hi=owner.time,t
            lower=None
            for _ in range(lim['event_bisections']):
                mid=.5*(lo+hi)
                b=branch(mid,released,'INTERIOR_BISECTION')
                valid=feasible_cone(b,released,owner.mu,lim,strict=True)
                if valid: lo,lower=mid,b
                else: hi=mid
                emit(dict(stage='BRACKET',lower_s=lo,upper_s=hi,strict_feasible=valid,
                          feasible_lower_seen=lower is not None))
                if hi-lo<=lim['event_time_tolerance_s']: break
            if lower is None or hi-lo>lim['event_time_tolerance_s']:
                raise CandidateFailure('onset interior event not bracketed',c)
            c=lower
            c.update(event_kind='CONE_LIMIT',event_bracket_s=[lo,hi])
        c.update(released=released.copy(),mode_trace=[eligible.tolist(),released.tolist()])
        owner.commit(c)
        emit(dict(stage='ACCEPTED',time_s=owner.time,step=owner.count,event_kind=c['event_kind']))
        if on_accept is not None: on_accept()
        if owner.time==t: return
    raise CandidateFailure('bounded onset events exhausted')


class OnsetBlock(RecoveryBlock):
    def __init__(self):
        super().__init__()
        cfg=json.loads(CONFIG.read_text(encoding='utf-8'))
        if cfg['variant']!='Q7_BOUNDARY_ELIGIBILITY_V1' or cfg['maximum_iterations']!=14 or cfg['production_authorized'] is not False:
            raise ValueError('undeclared onset variant')
        self.pfcfg['onset_variant']=cfg
        self.trace_sink=None

    def advance(self,t,on_accept=None):
        if any(self.time<k-1e-14 and t>k+1e-14 for k in (.05,.1,.15)):
            raise ValueError('requested interval crosses Q7 program corner')
        records=[]; before_time=self.time; before_step=self.count; error=None
        try:
            return advance_onset(self,t,on_accept,records.append)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            error=str(exc)
            raise
        finally:
            if self.trace_sink is not None:
                self.trace_sink(dict(classification='ATTEMPT_TRACE_NOT_RESTART',restart_allowed=False,
                    requested_time_s=t,base_time_s=before_time,base_step=before_step,
                    last_accepted_time_s=self.time,last_accepted_step=self.count,error=error,trace=records))


def main():
    if '--output' not in sys.argv:
        raise ValueError('explicit --output PATH required')
    output=Path(sys.argv[sys.argv.index('--output')+1])
    def construct():
        cfg=json.loads(CONFIG.read_text(encoding='utf-8'))
        provenance=dict(variant=cfg['variant'],production_authorized=False,configuration=cfg,
            source_sha256={s:parent.sha(ROOT/s) for s in EXTRA_SOURCES},
            baseline_source_sha256={s:parent.sha(ROOT/s) for s in parent.SOURCES})
        with (output/'onset_provenance.json').open('x',encoding='utf-8') as f:
            f.write(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
        block=OnsetBlock()
        def save_trace(data):
            name=f"attempt_{data['base_step']:04d}_{data['requested_time_s']:.17g}.json"
            with (output/name).open('x',encoding='utf-8') as f:
                f.write(json.dumps(data,indent=2,allow_nan=False)+'\n')
        block.trace_sink=save_trace
        return block
    with patch.object(parent,'PlasticFrictionBlock',construct):
        return parent.main()


if __name__=='__main__':
    raise SystemExit(main())
