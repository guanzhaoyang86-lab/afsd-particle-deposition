"""Integrate the verified predictor with frozen Onset events/atomic commit."""
import json
from pathlib import Path
import sys
from unittest.mock import patch

import probe_process_plastic_friction as parent
from probe_process_plastic_friction_onset import OnsetBlock, EXTRA_SOURCES as ONSET_SOURCES
from probe_process_feasible_predictor import solve_with_feasible_predictor

ROOT=parent.ROOT
CONFIG=ROOT/'configs/benchmarks/process_predictor_onset.json'
EXTRA_SOURCES=['scripts/probe_process_predictor_onset.py',
    'configs/benchmarks/process_predictor_onset.json',
    'scripts/probe_process_feasible_predictor.py']


def configuration():
    cfg=json.loads(CONFIG.read_text(encoding='utf-8'))
    if (cfg['variant']!='Q7_PREDICTOR_ONSET_V1' or cfg['production_authorized'] is not False
            or cfg['parent']!='configs/benchmarks/process_plastic_friction_onset.json'
            or type(cfg['maximum_total_corrections']) is not int
            or cfg['maximum_total_corrections']!=14
            or cfg['run_timeout_s']!=180 or cfg['permitted_cases']!=[4,8]):
        raise ValueError('undeclared predictor/onset integration')
    return cfg


class PredictorOnsetBlock(OnsetBlock):
    def __init__(self):
        super().__init__()
        self.pfcfg['predictor_variant']=configuration()

    def branch(self,t,released,maximum_iterations=None):
        maximum=14 if maximum_iterations is None else maximum_iterations
        return solve_with_feasible_predictor(self,t,released,maximum_iterations=maximum)


def main():
    if '--output' not in sys.argv:
        raise ValueError('explicit --output PATH required')
    output=Path(sys.argv[sys.argv.index('--output')+1])
    def construct():
        cfg=configuration()
        provenance=dict(variant=cfg['variant'],production_authorized=False,configuration=cfg,
            source_sha256={s:parent.sha(ROOT/s) for s in EXTRA_SOURCES},
            onset_source_sha256={s:parent.sha(ROOT/s) for s in ONSET_SOURCES},
            baseline_source_sha256={s:parent.sha(ROOT/s) for s in parent.SOURCES},
            note='PredictorOnsetBlock overrides branch ONLY; unchanged Onset advance and Q7 commit.')
        with (output/'predictor_provenance.json').open('x',encoding='utf-8') as f:
            f.write(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
        block=PredictorOnsetBlock()
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
