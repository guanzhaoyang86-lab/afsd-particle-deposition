"""Render a new figure with the upper-span caption at the left of panel b.

The established plotting script and all data remain unchanged. This adapter
only relocates one explanatory text object to keep the final orange data point
visible as the plotted values grow. It never imports or runs a physics solver.
"""
from pathlib import Path
from unittest.mock import patch
import importlib.util
import json
from matplotlib.axes import Axes

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('qualified_length_plot', BASE/'plot_comparison.py')
plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plot)
output = BASE/'plots_final_11p5625ms_layout_v2'
plot.require(not output.exists(), 'new layout output directory required')
records = plot.source_comparisons()
plot.require(len(records) == 10 and abs(records[-1]['data']['time_s']-.0115625)<1e-13,
             'this layout receipt is for the final 10-point 11.5625 ms comparison')
original_text = Axes.text


def relocated_text(axis, x, y, text, *args, **kwargs):
    if isinstance(text, str) and text.startswith('Different free spans (initial z):'):
        x = .02
        kwargs['ha'] = 'left'
    return original_text(axis, x, y, text, *args, **kwargs)


output.mkdir()
with patch.object(Axes, 'text', relocated_text):
    summary = plot.make_plot(records, output)
summary['layout_adapter_sha256'] = plot.sha(__file__)
summary['layout_only_change'] = 'Panel b span caption moved from upper right to upper left; no data, scales, or underlying plotting source changed.'
(output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
print(json.dumps(dict(output=str(output), matched_time_points=len(records),
                     latest_time_ms=summary['latest_time_s']*1e3), indent=2))
