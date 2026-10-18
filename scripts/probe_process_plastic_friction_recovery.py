"""Q7R: one explicit 14-correction budget, unchanged physics/event/commit.

Use only as a single-process companion. The scoped constructor injection reuses
the frozen Q7 CLI/output implementation; no parent file is edited.
"""
import copy
import json
from pathlib import Path
from unittest.mock import patch

import probe_process_plastic_friction as parent

ROOT = parent.ROOT
CONFIG = ROOT/'configs/benchmarks/process_plastic_friction_q7_recovery.json'
EXTRA_SOURCES = ['scripts/probe_process_plastic_friction_recovery.py',
                 'configs/benchmarks/process_plastic_friction_q7_recovery.json']


def recovery_config():
    cfg = json.loads(CONFIG.read_text(encoding='utf-8'))
    if (cfg['variant'] != 'Q7R_BOUNDED_NEWTON14' or cfg['production_authorized'] is not False
            or cfg['maximum_iterations']['value'] != 14 or cfg['maximum_iterations']['unit'] != 'count'
            or cfg['baseline_maximum_iterations'] != 12 or cfg['run_timeout_s'] != 60):
        raise ValueError('undeclared recovery budget/scope')
    return cfg


class RecoveryBlock(parent.PlasticFrictionBlock):
    def __init__(self):
        super().__init__()
        cfg = recovery_config()
        if self.qcfg['numerical']['maximum_iterations'] != cfg['baseline_maximum_iterations']:
            raise ValueError('parent iteration budget changed')
        self.qcfg = copy.deepcopy(self.qcfg)
        self.qcfg['numerical']['maximum_iterations'] = cfg['maximum_iterations']['value']
        self.pfcfg = copy.deepcopy(self.pfcfg)
        self.pfcfg['recovery_variant'] = cfg

    def branch(self, t, released, maximum_iterations=None):
        maximum = 14 if maximum_iterations is None else maximum_iterations
        if type(maximum) is not int or not 0 <= maximum <= 14:
            raise ValueError('Q7R requires integer correction budget in [0,14]')
        return super().branch(t, released, maximum_iterations=maximum)


def main():
    # Parent argparse remains the only CLI parser and creates a fresh directory.
    # During construction that directory and the original provenance exist.
    import sys
    args = sys.argv[1:]
    if '--output' not in args:
        raise ValueError('explicit --output PATH is required')
    output = Path(args[args.index('--output')+1])
    cfg = recovery_config()

    def construct():
        sidecar = dict(variant=cfg['variant'], production_authorized=False,
            source_sha256={s:parent.sha(ROOT/s) for s in EXTRA_SOURCES},
            baseline_source_sha256={s:parent.sha(ROOT/s) for s in parent.SOURCES},
            configuration=cfg,
            note='Supplement to original Q7 source inventory; actual constructor RecoveryBlock, budget14; unchanged event and commit.')
        with (output/'recovery_provenance.json').open('x', encoding='utf-8') as f:
            f.write(json.dumps(sidecar, indent=2, allow_nan=False)+'\n')
        return RecoveryBlock()

    # Explicit and temporary dependency injection; not a persistent backend swap.
    with patch.object(parent, 'PlasticFrictionBlock', construct):
        return parent.main()


if __name__ == '__main__':
    raise SystemExit(main())
