"""Real solved fields only; upstream NPZ checkpoints, project VTK time series."""
import hashlib
import json
from pathlib import Path
import numpy as np
from pysph.solver.utils import dump, load
from afsd_scaffold.vtk_legacy import write_point_cloud
from .mechanics import positions, velocities, pressure, volumetric_energy
from .material import equivalent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checkpoint(path, pa, state):
    # A hard wall-clock stop must not truncate the previous accepted checkpoint.
    target = Path(path)
    pending = target.with_name(target.stem+'.pending'+target.suffix)
    dump(str(pending), [pa], state, detailed_output=True, compress=True)
    pending.replace(target)


def restore(path):
    data = load(str(path))
    state = data['solver_data']
    if state.get('classification') == 'REJECTED' or state.get('restart_allowed') is False:
        raise ValueError('REJECTED candidate cannot be restored as an accepted physical state')
    return data['arrays']['solid'], data['solver_data']


def vtk(path, pa, p, bulk, time, classification='DEBUG'):
    S = pa.S.reshape(-1, 3, 3)
    scalar = {key: getattr(pa, key).tolist() for key in ('particle_id', 'material_id', 'component_id', 'tracer_type', 'wholebar_tag', 'actuator')}
    scalar.update(temperature_K=pa.temperature.tolist(), density_kg_per_m3=pa.rho.tolist(), pressure_Pa=pressure(pa.rho, p['rho0'], bulk).tolist(), equivalent_stress_Pa=equivalent(S).tolist(), equivalent_plastic_strain=pa.ep.tolist(), equivalent_deviatoric_strain_rate_per_s=pa.strain_rate.tolist(), equivalent_plastic_strain_rate_per_s=pa.ep_rate.tolist(), mass_kg=pa.m.tolist(), smoothing_length_m=pa.h.tolist())
    scalar['internal_energy_J_per_kg'] = (p['cp']*(pa.temperature-p['room'])+volumetric_energy(pa.rho, p['rho0'], bulk)+np.sum(S*S, axis=(1, 2))/(4*p['shear']*pa.rho)).tolist()
    sigma = S-pressure(pa.rho, p['rho0'], bulk)[:, None, None]*np.eye(3)
    for a, b, name in ((0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'), (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')):
        scalar['stress_'+name+'_Pa'] = sigma[:, a, b].tolist()
        scalar['strain_rate_'+name+'_per_s'] = pa.D.reshape(-1, 3, 3)[:, a, b].tolist()
    vectors = dict(velocity_m_per_s=velocities(pa).tolist(), initial_position_m=np.column_stack((pa.x0, pa.y0, pa.z0)).tolist())
    if time == 0:
        # Derivatives have not been evaluated at initialization. Do not publish
        # scratch zero arrays as computed rates (the actuator already moves).
        scalar = {k: v for k, v in scalar.items() if 'rate' not in k}
    write_point_cloud(path, positions(pa).tolist(), scalar, vectors, title=f'{classification}_THERMOMECHANICAL_SPH time_s={time:.17g}; NOT_PAPER_VALIDATION')
    return list(scalar)+list(vectors)+['position_m']
