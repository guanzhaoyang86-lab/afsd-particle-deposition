"""Independent hot-domain factories and geometry-only checks, never acceptance.

No old checkpoint is loaded, no source is patched, and no process controller
or old spatial-transition function is used to invent an independent parent.
The CLI does not call NNPS, residual, material, thermal or contact-force maps.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import math
import subprocess
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:
    sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import discrete_process
from afsd_backend.pysph.quasistatic_process import pose
from afsd_backend.pysph.quasistatic_contact import QuasistaticContact
from afsd_backend.pysph.output import checkpoint
from pysph.solver.utils import load
from afsd_backend.pysph.interval_acceptance import (
    validate_interval_energy_state, validate_bore_crossing_state)
from batched_contact import BatchedChainRuleContact
from visible_quasistatic import VisibleSupportTrial, CLASSIFICATION as VISIBLE_OPERATOR_CLASS

CLASSIFICATION = 'EXPANDED_HOT_QS_VISIBLE_SUPPORT_APPROXIMATE_V1'
SCOPE = 'INDEPENDENT_EXPANDED_HOT_450C_DOMAIN_INITIALIZATION_V1'
BACKEND_COMMIT = '1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5'
PARAMETER_SOURCE = ROOT/'outputs/contact_material_search_20261004/root/resolved_config.json'
PARAMETER_SHA = '37a153e2df917070ed707a93abcda8709cc8468f82393fc6cc90e07071bbdfd7'
CONTRACT = ROOT/'docs/EXPANDED_HOT_DOMAIN_INITIALIZATION_CONTRACT_20261004.md'
WIDTHS = (.006, .008)


def require(value, message):
    if not bool(value): raise ValueError(message)


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k:clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    return value


def write_new(path, value):
    encoded = json.dumps(clean(value), indent=2, allow_nan=False)+'\n'
    with Path(path).open('x', encoding='utf-8') as stream: stream.write(encoded)


def backend_provenance():
    import pysph
    source = Path(pysph.__file__).resolve().parents[1]
    expected = (ROOT/'third_party/PySPH').resolve()
    require(source == expected, 'initializer requires the pinned local PySPH checkout')
    # Trust this one already resolved project checkout for this read-only git
    # invocation; do not change the user's global Git configuration.
    commit = subprocess.check_output(['git', '-c', 'safe.directory='+source.as_posix(),
        '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    require(commit == BACKEND_COMMIT, 'PySPH backend commit differs')
    return dict(name='PySPH', commit=commit, source_directory=source.relative_to(ROOT).as_posix())


def parameter_units():
    groups = {
        'm': 'dx width bar_length substrate_width substrate_depth bore_width tool_radius tool_top layer marker_low marker_high actuator_depth contact_reference_spacing maximum_displacement_m',
        'K': 'temperature base_temperature tool_temperature room melt',
        'kg/m^3': 'rho0', 'Pa': 'shear a1 a2', 'W/(m K)': 'conductivity',
        'J/(kg K)': 'cp', 'W/(m^2 K)': 'contact_conductance', 'N/m': 'normal_stiffness',
        's^-1': 'reference_rate rate_floor', 'm/s': 'traverse feed lift_speed max_speed',
        'rad/s': 'omega', 's': 'wall_seconds dwell_duration release_time',
        'J': 'energy_absolute_limit', 'string': 'spatial_operator',
        '1': 'h_ratio poisson friction penalty_factor heat_partition plastic_heat_fraction a3 b1 b2 c1 c2 cfl density_limit penetration_limit energy_relative_limit',
    }
    return {key:unit for unit, keys in groups.items() for key in keys.split()}


def parameters(substrate_width_m):
    width = float(substrate_width_m)
    require(width in WIDTHS, 'only independent 6 mm and 8 mm domains are registered')
    require(sha(PARAMETER_SOURCE) == PARAMETER_SHA, 'immutable parameter source changed')
    base = json.loads(PARAMETER_SOURCE.read_text(encoding='utf-8'))['parameters']
    p = dict(base, substrate_width=width, bore_width=.002504)
    units = parameter_units()
    require(set(units) == set(p), 'every resolved parameter needs an explicit unit')
    entries = {key:dict(value=value, unit=units[key], source=PARAMETER_SOURCE.relative_to(ROOT).as_posix(),
        source_sha256=PARAMETER_SHA, status='INHERITED_ENGINEERING_INPUT_NOT_NEW_CALIBRATION') for key,value in p.items()}
    for key, reason in (('substrate_width', 'NEW_INDEPENDENT_DOMAIN_WIDTH_NOT_OLD_CHECKPOINT_EXTENSION'),
                        ('bore_width', 'INHERITED_UNLOADED_2UM_SIDE_GAP_ENGINEERING_APPROXIMATION')):
        entries[key].update(source=CONTRACT.relative_to(ROOT).as_posix(), source_sha256=sha(CONTRACT),
            original_value=base[key], status=reason)
    for key in ('a3', 'b1', 'b2', 'c1', 'c2', 'rate_floor'):
        entries[key]['limitation'] = 'Inherited FKS interpretation and low-rate assumptions are not paper-validated.'
    return p, dict(scope=SCOPE, parameter_source=PARAMETER_SOURCE.relative_to(ROOT).as_posix(),
        source_sha256=PARAMETER_SHA, contract_sha256=sha(CONTRACT), entries=entries,
        changed_parameter_names=sorted(k for k in p if p[k] != base[k]),
        old_state_read=False, parent_checkpoint=None, backend=backend_provenance(),
        mass_scale=1., time_scale=1., modulus_scale=1.)


def source_hashes():
    files = list((ROOT/'src/afsd_backend/pysph').glob('*.py'))
    files += [Path(__file__).resolve(), CONTRACT, PARAMETER_SOURCE]
    files += [ROOT/'scripts'/name for name in ('visible_quasistatic.py', 'visible_support.py',
        'visible_support_numerics.py', 'batched_contact.py', 'run_discrete_contact_process.py',
        'run_hot_quasistatic_probe.py', 'audit_solid_occluded_neighbors.py')]
    return {path.relative_to(ROOT).as_posix():sha(path) for path in sorted(set(files))}


def independent_initial_state(p, provenance=None, *, variant='visible'):
    require(variant in ('visible', 'discrete'), 'registered operator variant required')
    expected, evidence = parameters(p['substrate_width'])
    require(p == expected and (provenance is None or provenance == evidence), 'unregistered parameter/provenance change')
    pa, state = discrete_process.independent_initial_state(p)
    state.update(classification=CLASSIFICATION if variant == 'visible' else discrete_process.CLASSIFICATION,
        model_scope=SCOPE, operator_variant=variant, parameters=copy.deepcopy(p), parameter_provenance=evidence,
        source_sha256=source_hashes(), independent_initialization=True, restart_allowed=False,
        accepted=False, crossing_mass_kg=0., independent_domain_initialization=dict(
            kind=SCOPE, substrate_width_m=p['substrate_width'], initialization_time_s=0., initialization_step=0,
            parent_checkpoint=None, parent_checkpoint_sha256=None, old_state_read=False,
            inherited_accepted_time_s=0., id_namespace=f"independent_{int(round(p['substrate_width']*1000))}mm",
            spatial_operator_classification=VISIBLE_OPERATOR_CLASS if variant == 'visible' else discrete_process.CLASSIFICATION,
            physical_qualification='NOT_RUN', initialization_check='PARAMETERS_AND_ARRAYS_ONLY',
            scratch_rate_fields='UNEVALUATED_INITIAL_ZEROS_NOT_COMPUTED_RATES'))
    require('spatial_model_transition' not in state and 'branch_origin_sha256' not in state,
            'independent initialization cannot fabricate a prior spatial transition')
    validate_interval_energy_state(state['interval_energy']); validate_bore_crossing_state(state['bore_observation'])
    return pa, state


class ExpandedVisibleTrial(VisibleSupportTrial):
    """Diagnostic adapter; unchanged visible operator with an independent origin.

    This constructor has no NNPS or physical evaluation. Inherited evaluate and
    tangent become available only to a separately supervised cost/root caller.
    The old visible accepted_copy transition gate must not be used for this state.
    """
    def __init__(self, pa, state, p):
        require(state.get('classification') == CLASSIFICATION and state.get('model_scope') == SCOPE and
                state.get('independent_domain_initialization', {}).get('kind') == SCOPE and
                'spatial_model_transition' not in state, 'explicit independent expanded-visible origin required')
        require(state['parameters'] == p, 'trial parameters differ from independent state')
        adapter = dict(state, classification=discrete_process.CLASSIFICATION)
        discrete_process.DiscreteContactTrial.__init__(self, pa, adapter, p)
        self.state = state
        self.contact = BatchedChainRuleContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
            tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2,
            normal_stiffness=p['normal_stiffness'], tangential_stiffness=p['normal_stiffness'],
            friction=p['friction'], base_z=-p['substrate_depth'])


def make_trial(pa, state):
    p = state['parameters']
    if state['operator_variant'] == 'visible': return ExpandedVisibleTrial(pa, state, p)
    require(state['model_scope'] == SCOPE, 'independent expanded discrete scope required')
    return discrete_process.DiscreteContactTrial(pa, state, p)


def initial_geometry_and_field_audit(pa, state):
    p = state['parameters']; n = len(pa.x); dx = p['dx']
    x = np.column_stack((pa.x, pa.y, pa.z)); rod = pa.component_id == 1; substrate = ~rod
    require(state['time_s'] == 0. and state['step'] == 0 and state['restart_allowed'] is False and
            state['accepted'] is False and not state['interval_energy']['records'], 'initial audit cannot consume an evolved state')
    for key, prop in pa.properties.items():
        require(np.isfinite(prop.get_npy_array()).all(), 'nonfinite initial field: '+key)
    require(np.array_equal(pa.particle_id, np.arange(1,n+1)) and np.all(pa.material_id == 1), 'initial ID/material assignment differs')
    require(np.array_equal(x, np.column_stack((pa.x0,pa.y0,pa.z0))), 'initial position identity differs')
    for key in ('u','v','w','S','ep','plastic_work','return_dissipation','feed_crossed'):
        require(not np.any(getattr(pa,key)), 'initial unloaded/rest field differs: '+key)
    require(np.all(pa.temperature == 723.15) and np.all(pa.rho == p['rho0']) and
            np.all(pa.m == p['rho0']*dx**3) and np.all(pa.h == dx*p['h_ratio']), 'initial physical properties differ')
    bottom = substrate & (pa.z < -p['substrate_depth']+dx)
    top = rod & (pa.z > p['bar_length']-p['actuator_depth'])
    free = np.ones((n,3),bool); free[bottom] = False; free[top,2] = False
    require(np.array_equal(free,state['free']) and np.array_equal(top,pa.actuator.astype(bool)), 'prescribed axes differ')
    require(state['contact_history'].shape == (n,8,3) and not np.any(state['contact_history']), 'initial contact history differs')
    require(int(rod.sum()) == 150 and int(top.sum()) == 25, 'bar or driven top count changed')
    expected_substrate = int(round(p['substrate_width']/dx))**2*int(round(p['substrate_depth']/dx))
    require(int(substrate.sum()) == expected_substrate and len(np.unique(x,axis=0)) == n, 'substrate count or duplicate positions differ')
    minimum_squared = min(float(np.min(np.sum((x[i+1:]-x[i])**2,axis=1))) for i in range(n-1))
    geometry_scale = max(p['substrate_width'],p['bar_length'],p['tool_radius'])
    roundoff = 64*np.finfo(float).eps*geometry_scale
    require(math.sqrt(minimum_squared) >= dx-roundoff, 'initial nominal particle cells overlap')
    tool = QuasistaticContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'],particle_radius=dx/2,normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'],friction=p['friction'],base_z=-p['substrate_depth'])
    initial_pose = pose(p,0.); local=(x-initial_pose['origin'])@initial_pose['rotation']
    # _patches returns distances/projections only. Do not call geometry(),
    # which also calculates forces, or any contact/residual evaluate method.
    gaps=np.column_stack((np.array([tool._patches(point)[0] for point in local]), x[:,2]+p['substrate_depth']))
    signed_clearance=gaps-dx/2
    require(np.isfinite(gaps).all() and float(signed_clearance.min()) >= -roundoff,
            'initial particle shell enters rigid finite tool/base beyond roundoff')
    overlap=np.maximum(-signed_clearance,0.)
    initial_normal_potential=math.fsum((.5*p['normal_stiffness']*overlap**2).ravel())
    eligible=pa.particle_id[rod & (pa.z>p['layer'])].astype(int).tolist()
    require(state['bore_observation']['eligible_ids'] == eligible and len(eligible)==75 and
            not state['crossing_ids'] and state['crossing_mass_kg']==0., 'initial throughput eligibility differs')
    mass=math.fsum(pa.m); volume=math.fsum(pa.m/pa.rho); nf=int(free.sum())
    dense=dict(gradient_map_D_bytes=8*(9*n)*(3*n),full_force_tangent_J_bytes=8*(3*n)**2,
        free_matrix_bytes=8*nf**2,free_inverse_bytes=8*nf**2,thermal_matrix_bytes=8*n*n)
    return dict(status='PASS_INITIAL_ARRAYS_AND_GEOMETRY_ONLY_NOT_PHYSICAL_QUALIFICATION',
        classification=state['classification'],particle_count=n,substrate_count=int(substrate.sum()),rod_count=int(rod.sum()),
        bottom_fixed_particles=int(bottom.sum()),top_z_prescribed_particles=int(top.sum()),free_displacement_dof=nf,
        mass_kg=mass,volume_m3=volume,eligible_initial_rod_mass_kg=math.fsum(pa.m[np.isin(pa.particle_id,eligible)]),
        initially_lower_rod_mass_kg=math.fsum(pa.m[rod & (pa.z<=p['layer'])]),eligible_ids=eligible,
        initial_energy=dict(kinetic_J=0.,deviatoric_elastic_J=0.,volumetric_elastic_J=0.,plastic_work_J=0.,
            tangential_contact_potential_J=0.,normal_contact_potential_J=initial_normal_potential,
            sensible_heat_reference_K=p['room'],sensible_heat_J=math.fsum(pa.m*p['cp']*(pa.temperature-p['room'])),
            force_equilibrium='NOT_EVALUATED'),
        geometry=dict(minimum_particle_centre_distance_m=math.sqrt(minimum_squared),
            minimum_signed_shell_clearance_m=float(signed_clearance.min()),roundoff_bound_m=roundoff,
            roundoff_shell_overlap_max_m=float(overlap.max()),base_z_m=-p['substrate_depth'],tool_pose=initial_pose,
            feature_names=list(tool.features),square_bore_open=True,tool_height_m=p['tool_top']-p['layer'],
            initial_shoulder_edge_margin_m=p['substrate_width']/2-p['tool_radius'],kernel_support_radius_m=2*dx*p['h_ratio'],
            nominal_unbuffered_static_plate_coverage_time_s=(p['substrate_width']/2-p['tool_radius'])/p['traverse']),
        theoretical_float64_dense_allocations=dense,allocated_dense_matrices=False,
        physical_solver_evaluations=0,nnps_evaluations=0,contact_force_evaluations=0,
        accepted_steps=0,accepted_time_advanced_s=0.,spatial_domain_sensitivity='NOT_RUN',
        uncomputed_fields=['D','strain_rate','ep_rate','ax','ay','az','au','av','aw','arho','ae','e'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width-mm',type=int,choices=(6,8),required=True)
    parser.add_argument('--variant',choices=('visible','discrete'),default='visible')
    parser.add_argument('--output',required=True)
    args=parser.parse_args(); output=Path(args.output)
    output=(ROOT/output).resolve() if not output.is_absolute() else output.resolve()
    require(output.is_relative_to(ROOT) and not output.exists(), 'new output inside project required')
    output.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    p,provenance=parameters(args.width_mm*.001);pa,state=independent_initial_state(p,provenance,variant=args.variant)
    audit=initial_geometry_and_field_audit(pa,state)
    # Diagnostic constructor and prescribed axes only; evaluate/tangent unused.
    trial=make_trial(pa,state); prescribed=trial.prescribed(1.)
    require(np.array_equal(prescribed[:,2],-p['feed']*pa.actuator) and not np.any(prescribed[:,:2]),
            'lightweight constructor changed bar-top feed axes')
    state['initial_array_geometry_audit']=audit
    snapshot=output/'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'; checkpoint(snapshot,pa,state)
    restored=load(str(snapshot)); read_pa=restored['arrays']['solid']; read_state=restored['solver_data']
    require(set(read_pa.properties)==set(pa.properties) and clean(read_state)==clean(state),
            'independent snapshot field schema/state roundtrip differs')
    for key,prop in pa.properties.items():
        require(np.array_equal(prop.get_npy_array(),read_pa.properties[key].get_npy_array()),
                'independent snapshot particle field roundtrip differs: '+key)
    write_new(output/'resolved_parameters.json',dict(parameters=p,provenance=provenance))
    write_new(output/'initialization.json',state)
    fields={key:prop.get_npy_array().copy() for key,prop in pa.properties.items()}
    np.savez_compressed(output/'initial_particle_fields.npz',accepted=False,restart_allowed=False,
        time_s=0.,step=0,free_axes=state['free'],contact_history=state['contact_history'],**fields)
    for path,digest in state['source_sha256'].items(): require(sha(ROOT/path)==digest,'source changed during initialization')
    report=dict(audit,command=[sys.executable,*sys.argv],exit_code=0,wall_seconds=time.perf_counter()-begin,
        complete_particle_and_state_checkpoint_roundtrip='PASS',
        initialization_checkpoint=snapshot.relative_to(ROOT).as_posix(),initialization_sha256=sha(snapshot),
        source_sha256=state['source_sha256'],output_sha256={path.name:sha(path) for path in sorted(output.iterdir()) if path.is_file()})
    write_new(output/'summary.json',report)
    print(json.dumps(clean({key:report[key] for key in ('status','particle_count','free_displacement_dof',
        'mass_kg','theoretical_float64_dense_allocations','initialization_sha256','wall_seconds')}),indent=2))
    return 0


if __name__ == '__main__': raise SystemExit(main())
