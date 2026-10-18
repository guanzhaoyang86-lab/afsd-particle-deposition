"""Explicit affordable horizon adapter; same coupled physics, old cost untouched."""
from pathlib import Path
import argparse
import json
import math
import numpy as np
import run_coupled_afsd as base

ROOT=base.ROOT
ORIGINAL_SETTINGS=base.settings
ORIGINAL_SOURCES=list(base.SOURCES)
EXTRA=['scripts/run_coupled_short_window.py','docs/COUPLED_AFSD_SHORT_WINDOW.md']
PROBE=ROOT/'outputs/coupled_afsd_20261002/cost_probe/manifest.json'


def original_cost():
    m=json.loads(PROBE.read_text(encoding='utf-8'))
    cfg,p,resolved=ORIGINAL_SETTINGS()
    if m['status']!='COMPLETE_COST_PROBE' or m['accepted_steps']!=5:
        raise ValueError('original accepted cost probe required')
    if m['source_sha256']!={s:base.sha(ROOT/s) for s in ORIGINAL_SOURCES}:
        raise ValueError('original runtime source changed')
    if m['configuration_sha256']!=base.fingerprint(resolved):raise ValueError('original configuration changed')
    for name,h in m['artifact_sha256'].items():
        if Path(name).name!=name or base.sha(PROBE.parent/name)!=h:raise ValueError('original cost artifact changed')
    return m,cfg,p,resolved


def affordable_choice(m,cfg):
    dt=min(m['accepted_dt_s']);rate=max(m['step_wall_seconds'])
    if not np.isfinite([dt,rate]).all() or dt<=0 or rate<=0:raise ValueError('invalid cost measurement')
    unit=1e-7;capacity=(cfg['forecast_launch_limit_seconds']-cfg['forecast_overhead_seconds'])/(cfg['forecast_multiplier']*rate)
    count=math.floor(capacity*dt/unit)
    while count>0:
        end=count*unit;steps=math.ceil(end/dt)
        predicted=cfg['forecast_multiplier']*rate*steps+cfg['forecast_overhead_seconds']
        if predicted<=cfg['forecast_launch_limit_seconds']:
            return dict(end_time_s=end,estimated_steps=steps,forecast_wall_seconds=predicted,permitted=True)
        count-=1
    raise ValueError('not even one 100ns quantum affordable; NOT_RUN')


def short_settings():
    m,cfg,p,resolved=original_cost();choice=affordable_choice(m,cfg)
    # Only schedule metadata changes. p (all physical inputs) is untouched.
    new_cfg=dict(cfg,registered_end_times_s=[choice['end_time_s']])
    result=dict(resolved,explicit_configuration=new_cfg,
        original_explicit_configuration=cfg,
        short_window_budget_override=dict(classification='AFFORDABLE_SHORT_WINDOW_NOT_ORIGINAL_2US_COMPLETION',
            contract='docs/COUPLED_AFSD_SHORT_WINDOW.md',choice=choice,
            origin_probe_manifest_sha256=base.sha(PROBE),physical_parameters_unchanged=True))
    return new_cfg,p,result


def short_budget_report(probe_manifest):
    if Path(probe_manifest).resolve()!=PROBE.resolve():raise ValueError('only registered cost evidence permitted')
    m,cfg,_,_=original_cost();new_cfg,_,resolved=short_settings();choice=affordable_choice(m,cfg)
    return dict(status='READY_FOR_ONE_REGISTERED_TRIAL',production_authorized=False,
        probe_manifest_path=PROBE.relative_to(ROOT).as_posix(),probe_manifest_sha256=base.sha(PROBE),
        source_sha256={s:base.sha(ROOT/s) for s in ORIGINAL_SOURCES+EXTRA},
        configuration_sha256=base.fingerprint(resolved),choice=choice,
        evidence_status='TRANSFER_ORIGINAL_MEASURED_COST_SAME_PHYSICAL_PARAMETERS_NOT_A_NEW_PROBE',
        physical_parameters_unchanged=True,original_long_horizons_status='NOT_RUN_COST_EXCEEDS_BUDGET',
        internal_limit_s=cfg['internal_wall_seconds'],external_limit_s=cfg['external_wall_seconds'],
        method='largest 100ns quantum satisfying unchanged 1.5*slowest_step*ceil(end/dt)+30<=450s')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    g=p.add_mutually_exclusive_group(required=True);g.add_argument('--select-budget',action='store_true');g.add_argument('--budget',type=Path)
    a=p.parse_args()
    if a.select_budget:
        if a.output.exists():raise FileExistsError('fresh short budget output required')
        r=short_budget_report(PROBE);base.write_json(a.output,r);print(json.dumps(r,allow_nan=False));return 0
    # Explicit in-process entry-point binding, not a change to solver/physics.
    base.settings=short_settings;base.SOURCES=ORIGINAL_SOURCES+EXTRA;base.budget_report=short_budget_report
    return base.run(a.output,budget=a.budget)


if __name__=='__main__':raise SystemExit(main())
