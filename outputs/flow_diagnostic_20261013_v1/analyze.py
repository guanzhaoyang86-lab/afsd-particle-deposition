"""Passive, hash-checked flow diagnostics of complete saved endpoints."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
import run_gpu_nonlinear_resume_20261011 as runner
from afsd_backend.pysph.mechanics import bulk_modulus, pressure

start = runner.BASE / 'event_resume_8h_20261011_v1/attempt_0039'
campaign = runner.BASE / 'nonlinear_resume_20261011_v1'
stages = [start] + sorted(campaign.glob('attempt_[0-9][0-9][0-9][0-9]'))
records, sources = [], []
first = None
for stage in stages:
    manifest_path = stage / 'manifest.json'
    if not manifest_path.exists():
        continue
    manifest = runner.read(manifest_path)
    if 'screened_continuation_checkpoint' not in manifest:
        continue
    path = ROOT / manifest['screened_continuation_checkpoint']
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == manifest['screened_continuation_sha256']
    pa, state = runner.load(path)
    ids = np.asarray(pa.particle_id).copy()
    x = np.column_stack((pa.x, pa.y, pa.z))
    velocity = np.column_stack((pa.u, pa.v, pa.w))
    bar = pa.component_id == 1
    if first is None:
        first = dict(ids=ids, x=x.copy(), m=pa.m.copy(), z0=pa.z0.copy(), bar=bar.copy(),
            ep=pa.ep.copy(), plastic_work=pa.plastic_work.copy())
        layers = np.unique(pa.z0[bar])
        tracked = list(map(int, state['bore_observation']['crossing_ids']))
    assert np.array_equal(ids, first['ids']) and np.array_equal(pa.m, first['m'])
    assert np.array_equal(pa.z0, first['z0']) and np.array_equal(bar, first['bar'])
    if records:
        assert state['time_s'] > records[-1]['time_s']
    params = state['parameters']
    pose = state['machine_history'][-1]
    assert np.allclose(pose['origin1_m'], [0., 0., .0015], atol=1e-14, rtol=0)
    hydro = pressure(pa.rho, params['rho0'], bulk_modulus(params['shear'], params['poisson']))
    sigma = pa.S.reshape(-1, 3, 3) - hydro[:, None, None] * np.eye(3)
    layer_rows = []
    for z0 in layers:
        mask = bar & (pa.z0 == z0)
        xy = x[mask, :2]
        centered = xy - xy.mean(axis=0)
        xy0 = first['x'][mask, :2]
        centered0 = xy0 - xy0.mean(axis=0)
        radial = np.linalg.norm(xy, axis=1)
        radial_velocity = np.sum(xy * velocity[mask, :2], axis=1) / np.maximum(radial, 1e-30)
        layer_rows.append(dict(
            initial_z_mm=float(z0*1e3), count=int(mask.sum()),
            current_mean_z_mm=float(x[mask, 2].mean()*1e3),
            mean_dz_um=float((x[mask, 2]-first['x'][mask, 2]).mean()*1e6),
            mean_w_mm_s=float(velocity[mask, 2].mean()*1e3),
            mean_radial_velocity_mm_s=float(radial_velocity.mean()*1e3),
            centered_rms_radius_mm=float(np.sqrt(np.sum(centered**2, axis=1).mean())*1e3),
            centered_rms_radius_change_um=float((np.sqrt(np.sum(centered**2, axis=1).mean())-
                np.sqrt(np.sum(centered0**2, axis=1).mean()))*1e6),
            axis_offset_mm=float(np.linalg.norm(xy.mean(axis=0))*1e3),
            mean_pressure_MPa=float(hydro[mask].mean()*1e-6),
            mean_sigma_zz_MPa=float(sigma[mask, 2, 2].mean()*1e-6),
            mean_ep=float(pa.ep[mask].mean()),
            mean_ep_increment=float((pa.ep[mask]-first['ep'][mask]).mean()),
            temperature_mean_K=float(pa.temperature[mask].mean()),
            active_plastic_count_last_step=int(np.count_nonzero(pa.ep_rate[mask] > 0.)),
            plastic_work_increment_J=float((pa.plastic_work[mask]-first['plastic_work'][mask]).sum()),
        ))
    cross_mask = np.isin(ids, state['bore_observation']['crossing_ids'])
    local = (x - np.array(pose['origin1_m'])) @ np.array(pose['rotation1'])
    under_geometry = (x[:, 2] > 0.) & (local[:, 2] < 0.) & (
        np.max(np.abs(local[:, :2]), axis=1) >= params['bore_width']/2) & (
        np.linalg.norm(x[:, :2], axis=1) <= params['tool_radius'])
    under = cross_mask & under_geometry
    tracks = []
    for pid in tracked:
        k = int(np.flatnonzero(ids == pid)[0])
        tracks.append(dict(particle_id=pid, xyz_mm=(x[k]*1e3).tolist(),
            radius_mm=float(np.linalg.norm(x[k,:2])*1e3),
            velocity_mm_s=(velocity[k]*1e3).tolist(),
            pressure_MPa=float(hydro[k]*1e-6), sigma_zz_MPa=float(sigma[k,2,2]*1e-6)))
    regions = {}
    for name, mask in [('original_upper_free_bar',bar & (pa.z0 > .0025) & (pa.actuator==0)),
                       ('original_bore_bar',bar & (pa.z0 > .0015) & (pa.z0 < .0025)),
                       ('original_lower_bar',bar & (pa.z0 < .0015)),
                       ('actuator',pa.actuator!=0),('substrate',~bar)]:
        regions[name] = dict(count=int(mask.sum()),
            temperature_min_K=float(pa.temperature[mask].min()),
            temperature_max_K=float(pa.temperature[mask].max()),
            active_plastic_count_last_step=int(np.count_nonzero(pa.ep_rate[mask]>0.)),
            mean_ep_increment=float((pa.ep[mask]-first['ep'][mask]).mean()),
            plastic_work_increment_J=float((pa.plastic_work[mask]-first['plastic_work'][mask]).sum()))
    records.append(dict(time_s=state['time_s'], step=int(state['step']), regions=regions,
        stage=str(stage.relative_to(ROOT)), layers=layer_rows, tracked_crossed=tracks,
        crossed_count=int(cross_mask.sum()), under_shoulder_count=int(under.sum()),
        all_bar_under_shoulder_count=int(np.count_nonzero(bar & under_geometry)),
        original_lower_bar_under_shoulder_count=int(np.count_nonzero(bar & (pa.z0<.0015) & under_geometry))))
    sources.append(dict(checkpoint=str(path.relative_to(ROOT)), sha256=digest))

assert len(records) == 84
result = dict(classification='PASSIVE_COMPLETE_SAVED_STATE_FLOW_DIAGNOSTIC',
    physical_steps_executed=0, sources=sources, records=records,
    conventions=dict(positive_z='upward', pressure='K*(rho/rho0-1); positive compression',
        sigma_zz='Szz-pressure; tension positive',
        layers='25 persistent particles grouped by their initial z0',
        velocity='saved displacement/dt of the last fine step; not a long-time average',
        radial_size='RMS about each layer centroid; not a cross-sectional area or mass flux'),
    limitations=['Saved endpoints do not resolve unsaved trajectories.',
        'Layer averages and stresses do not establish a causal model defect.',
        'These coarse engineering-approximate states do not establish spatial convergence.'])
(OUT/'diagnostic.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')

print("Readback complete:",len(records),"saved states; output:",OUT/"diagnostic.json")
