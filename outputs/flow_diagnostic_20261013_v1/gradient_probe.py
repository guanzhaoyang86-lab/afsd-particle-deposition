"""Read-only frozen saved-graph probes, never physical imposed trajectories.

No solver, NNPS rebuild, constitutive update, checkpoint write or eigensolve.
Elastic energies below characterize a diagnostic operator with the saved
volumes; they are not actual material work or the prestressed plastic tangent.
"""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
import run_gpu_nonlinear_resume_20261011 as resume
from afsd_backend.pysph.quasistatic_process import _geometry
from afsd_backend.pysph.mechanics import bulk_modulus


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference(path):
    return dict(path=str(path.relative_to(ROOT)), sha256=sha(path))


def apply_sparse(u, volume, i, j, g, correction):
    """Sparse application of the saved gradient_map, no dense global matrix."""
    du = u[j] - u[i]
    gi = np.einsum('nji,nj->ni', correction[i], g)
    gj = np.einsum('nji,nj->ni', correction[j], g)
    result = np.zeros((len(u), 3, 3))
    np.add.at(result, i, volume[j, None, None] * du[:, :, None] * gi[:, None, :])
    np.add.at(result, j, volume[i, None, None] * du[:, :, None] * gj[:, None, :])
    return result


def metrics(u, H, x, V, i, j, mask, shear, bulk):
    sym = .5 * (H + H.swapaxes(1, 2))
    trace = np.trace(H, axis1=1, axis2=2)
    dev = sym - trace[:, None, None] * np.eye(3) / 3
    energy = V * (shear * np.sum(dev * dev, axis=(1, 2)) + .5 * bulk * trace**2)
    rms = lambda a: float(np.sqrt(np.sum(V[mask] * a[mask]) / V[mask].sum()))
    edges = mask[i] & mask[j]
    dr = x[j] - x[i]
    du = u[j] - u[i]
    weight = V[i] * V[j]
    edge_rms = float(np.sqrt(np.sum(weight[edges] * np.sum(du[edges]**2, axis=1) /
        np.sum(dr[edges]**2, axis=1)) / weight[edges].sum()))
    sym_rms = rms(np.sum(sym**2, axis=(1, 2)))
    result = dict(particles=int(mask.sum()), within_region_edges=int(edges.sum()),
        volume_m3=float(V[mask].sum()), displacement_volume_rms_m=rms(np.sum(u**2, axis=1)),
        neighbor_difference_rms_gradient=edge_rms,
        full_gradient_volume_rms=rms(np.sum(H**2, axis=(1, 2))),
        symmetric_gradient_volume_rms=sym_rms,
        axial_gradient_volume_rms=rms(H[:, 2, 2]**2),
        axial_gradient_volume_mean=float(np.sum(V[mask]*H[mask, 2, 2])/V[mask].sum()),
        diagnostic_unstressed_elastic_energy_J=float(energy[mask].sum()),
        sym_gradient_to_neighbor_difference_ratio=sym_rms/edge_rms if edge_rms else None,
        normalized_elastic_response=float(energy[mask].sum() /
            (V[mask].sum()*shear*edge_rms**2)) if edge_rms else None)
    return result


def alternating_projection(actual, alternating, x, V, mask):
    """Weighted descriptive axial fit, not a force/energy mode decomposition."""
    design = np.column_stack((np.ones(mask.sum()), x[mask]))
    weight = np.sqrt(V[mask] / V[mask].mean())
    target = actual[mask, 2]
    pattern = alternating[mask, 2]
    def residual(values):
        coeff = np.linalg.lstsq(design*weight[:, None], values*weight, rcond=None)[0]
        return values - design@coeff
    target_residual, pattern_residual = residual(target), residual(pattern)
    inner = lambda a, b: float(np.sum(V[mask]*a*b))
    norm = inner(pattern_residual, pattern_residual)
    coefficient = inner(target_residual, pattern_residual)/norm
    total = inner(target_residual, target_residual)
    return dict(affine_removed='Weighted least squares of [1,x,y,z] in current coordinates; diagnostic only.',
        fitted_pattern_multiplier=coefficient,
        fitted_nominal_alternating_amplitude_m=coefficient*float(np.max(np.abs(pattern))),
        actual_nonaffine_axial_volume_rms_m=float(np.sqrt(total/V[mask].sum())),
        pattern_fraction_of_nonaffine_axial_squared_displacement=coefficient**2*norm/total if total else None,
        meaning='Displacement projection only, not fraction of strain energy/plastic work or proof of an hourglass mode.')


def frame(source):
    checkpoint = ROOT/source['checkpoint']
    assert sha(checkpoint) == source['sha256']
    pa, state = resume.load(checkpoint)
    proposal = checkpoint.parent/'proposal'
    graph_path = proposal/'visible_graph.npz'
    solved_path = proposal/'SOLVED_UNCOMMITTED_candidate.npz'
    motion_path = proposal/'motion.json'
    source_references = [reference(path) for path in (checkpoint, graph_path, solved_path, motion_path)]
    with np.load(graph_path, allow_pickle=False) as data:
        graph = {k: data[k].copy() for k in data.files}
    with np.load(solved_path, allow_pickle=False) as data:
        solved = {k: data[k].copy() for k in data.files}
    motion = read(motion_path)
    x = np.column_stack((pa.x, pa.y, pa.z))
    actual = solved['displacement_m']
    V, i, j, g, B = (graph[k] for k in ('volume', 'i', 'j', 'g', 'correction'))
    assert np.array_equal(pa.particle_id, graph['particle_id'])
    assert np.array_equal(pa.particle_id, solved['particle_id'])
    assert np.array_equal(x, graph['position1']) and np.array_equal(x, solved['position_m'])
    np.testing.assert_allclose(x - graph['position0'], actual, rtol=1e-10, atol=2e-18)
    np.testing.assert_allclose(V, pa.m/pa.rho, rtol=2e-13, atol=0)
    p = state['parameters']; dx = p['dx']; shear = p['shear']; bulk = bulk_modulus(shear, p['poisson'])
    bar = np.asarray(pa.component_id) == 1
    groups = dict(all_material=np.ones(len(x), bool), all_bar=bar,
        original_upper_free_bar=bar & (pa.z0 > .0025) & (pa.actuator == 0),
        original_bore_bar=bar & (pa.z0 > .0015) & (pa.z0 < .0025),
        original_lower_bar=bar & (pa.z0 < .0015))
    # Amplitudes are arbitrary small linear diagnostic inputs, never commands.
    amplitude = 0.5e-9
    layer_float = np.asarray(pa.z0)/dx - .5
    layer = np.rint(layer_float).astype(int)
    np.testing.assert_allclose(layer_float, layer, rtol=0, atol=1e-12)
    alternating = np.zeros_like(x)
    alternating[:, 2] = amplitude*np.where(layer % 2 == 0, 1., -1.)
    affine = np.zeros_like(x)
    epsilon = 2*amplitude/dx
    affine[:, 2] = epsilon*(x[:, 2] - np.sum(V*x[:, 2])/V.sum())
    actual_axial = np.zeros_like(x); actual_axial[:, 2] = actual[:, 2]
    fields = dict(affine_axial_reference=affine, alternating_initial_z_layers=alternating,
        actual_last_fine_step=actual, actual_last_fine_step_axial_only=actual_axial)
    rows = {}; validation = {}; gradients = {}
    for name, u in fields.items():
        H = apply_sparse(u, V, i, j, g, B)
        # Independent implementation of raw pair sums / rebuilt 3x3 inverse.
        reference_H, reference_B = _geometry(x, u, pa.m, pa.rho, i, j, g)
        np.testing.assert_allclose(B, reference_B, rtol=5e-12, atol=5e-12)
        scale = max(float(np.max(np.abs(H))), np.finfo(float).tiny)
        relative = float(np.max(np.abs(H-reference_H))/scale)
        assert relative < 1e-10
        gradients[name] = H
        validation[name] = dict(sparse_vs_raw_rebuilt_gradient_relative_max_error=relative)
        rows[name] = {key: metrics(u, H, x, V, i, j, mask, shear, bulk) for key, mask in groups.items()}
    expected = np.zeros((3, 3)); expected[2, 2] = epsilon
    affine_error = float(np.max(np.abs(gradients['affine_axial_reference']-expected))/epsilon)
    assert affine_error < 1e-10
    analytic_energy = V.sum()*(2*shear/3 + bulk/2)*epsilon**2
    measured_energy = rows['affine_axial_reference']['all_material']['diagnostic_unstressed_elastic_energy_J']
    affine_energy_error = float(abs(measured_energy/analytic_energy-1))
    assert affine_energy_error < 1e-10
    actual_D = .5*(gradients['actual_last_fine_step']+gradients['actual_last_fine_step'].swapaxes(1, 2))/motion['dt_s']
    saved_D = np.asarray(pa.D).reshape(-1, 3, 3)
    actual_error = float(np.max(np.abs(actual_D-saved_D))/max(np.max(np.abs(saved_D)), np.finfo(float).tiny))
    assert actual_error < 1e-10
    comparison = {}
    for key in groups:
        alt = rows['alternating_initial_z_layers'][key]
        control = rows['affine_axial_reference'][key]
        factor = alt['neighbor_difference_rms_gradient']/control['neighbor_difference_rms_gradient']
        comparison[key] = dict(affine_displacement_multiplier_for_equal_neighbor_difference=factor,
            alternating_to_equal_neighbor_affine_symmetric_gradient_ratio=alt['symmetric_gradient_volume_rms']/
                (control['symmetric_gradient_volume_rms']*factor),
            alternating_to_equal_neighbor_affine_elastic_energy_ratio=alt['diagnostic_unstressed_elastic_energy_J']/
                (control['diagnostic_unstressed_elastic_energy_J']*factor**2))
    assert all(sha(ROOT/item['path']) == item['sha256'] for item in source_references)
    return dict(time_s=state['time_s'], step=state['step'], last_fine_dt_s=motion['dt_s'],
        sources=source_references, input_sources_unchanged_after_probe=True,
        graph_edges=len(i), parameters=dict(dx_m=dx, shear_Pa=shear, bulk_Pa=bulk),
        probe_amplitude_m=amplitude, affine_reference_gradient=epsilon,
        controls=dict(affine_independent_gradient_relative_max_error=affine_error,
            affine_independent_elastic_energy_relative_error=affine_energy_error,
            actual_gradient_vs_saved_D_relative_max_error=actual_error, field_identity_checks=validation),
        field_metrics=rows, equal_neighbor_difference_comparisons=comparison,
        actual_axial_projection={key: alternating_projection(actual, alternating, x, V, groups[key])
            for key in ('all_bar', 'original_upper_free_bar', 'original_lower_bar')})


def main():
    start = time.perf_counter()
    source_path = OUT/'diagnostic.json'
    diagnostic = read(source_path)
    selected = [diagnostic['sources'][0], diagnostic['sources'][-1]]
    rows = [frame(source) for source in selected]
    result = dict(classification='READ_ONLY_FROZEN_SAVED_GRAPH_GRADIENT_RESPONSE_DIAGNOSTIC',
        physical_steps_executed=0, checkpoints_modified=False, nonlinear_solves=0,
        arbitrary_perturbations_are_physical_trajectories=False,
        source=reference(source_path), script=reference(Path(__file__).resolve()),
        formulas=dict(gradient='H_i=sum_j V_j (u_j-u_i) outer (B_i.T g_ij), from archived visible graph/volume/B.',
            affine='u_z=epsilon*(current z-volume centroid z), epsilon=2*A/dx; independent exact H_zz=epsilon.',
            alternating='u_z=A*(-1)^integer_layer for every material point, integer_layer=round(z0/dx-0.5); other axes zero.',
            neighbor_difference='sqrt(sum_e Vi*Vj*||uj-ui||^2/||xj-xi||^2 / sum_e Vi*Vj); within-region edges only.',
            region_gradient='Volume RMS over region particles, evaluated with the full graph including cross-region support.',
            energy='sum_i Vi*(G*dev(sym H):dev(sym H)+K*tr(H)^2/2); unstressed constant-elastic reference, not actual plastic energy.',
            matched_comparison='Affine field scaled so the selected region has identical neighbor-difference RMS; ratios compare resulting symmetric-gradient RMS and elastic reference energy.'),
        records=rows, elapsed_s=time.perf_counter()-start,
        limitations=['No empirical response threshold or hourglass PASS/FAIL classification is assigned.',
            'Layer-alternating perturbation is a highest reference-grid-frequency diagnostic, not a smooth continuum solution or imposed boundary-compatible physical perturbation.',
            'Some manufactured perturbations move fixed substrate/cap degrees of freedom. They test the spatial gradient only, not admissible free mechanical modes.',
            'Neighbor-difference normalization is descriptive and includes all directions. Full rigid rotations have nonzero edge differences but zero symmetric gradient.',
            'Two deformed graphs and their final fine-step increments do not establish spatial convergence, history-integrated errors, or a unique cause of absent deposition.',
            'The diagnostic elastic reference excludes prestress, plasticity, thermal and contact tangent contributions; positive or weak response is not physical stability certification.',
            'Actual-step axial projection removes a best affine trend and reports displacement variance only; nonaffine physical flow can also have such a component.'])
    (OUT/'gradient_probe.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(elapsed_s=result['elapsed_s'], records=[dict(time_s=r['time_s'],
        controls=r['controls'], comparisons=r['equal_neighbor_difference_comparisons'],
        projection=r['actual_axial_projection']) for r in rows]), indent=2))


if __name__ == '__main__':
    main()
