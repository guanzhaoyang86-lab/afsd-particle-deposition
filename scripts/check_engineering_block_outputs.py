"""Read back E2 artifacts/IDs and impulse identity; not physical accuracy PASS."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from pysph.solver.utils import load
from afsd_backend.pysph.output import restore
from afsd_scaffold.vtk_legacy import read_point_cloud

directory = Path(sys.argv[1])
names = ('n5_cfl010', 'n5_cfl005', 'n5_cfl0025', 'n7_cfl005', 'n9_cfl005')
errors = []
profiles = []
for name in names:
    case = directory/name
    report = json.loads((case/'report.json').read_text(encoding='utf-8'))
    for file, digest in report['source_sha256'].items():
        assert hashlib.sha256((ROOT/file).read_bytes()).hexdigest() == digest, file
    for file, digest in report['artifact_sha256'].items():
        assert hashlib.sha256((case/file).read_bytes()).hexdigest() == digest, file
    with np.load(case/'history.npz', allow_pickle=False) as history:
        vtk = read_point_cloud(case/'final.vtk')
        pa = load(str(case/'checkpoint.npz'))['arrays']['solid']
        np.testing.assert_array_equal(history['particle_id'], pa.particle_id)
        np.testing.assert_array_equal(history['particle_id'], vtk['scalars']['particle_id'])
        np.testing.assert_array_equal(history['position_m'][-1], np.column_stack((pa.x, pa.y, pa.z)))
        np.testing.assert_array_equal(history['position_m'][-1], vtk['points'])
        np.testing.assert_array_equal(history['velocity_m_s'][-1], vtk['vectors']['velocity_m_s'])
        np.testing.assert_array_equal(history['stress_Pa'][-1, :, 0, 0], vtk['scalars']['stress_xx_Pa'])
        # Post-hoc error localization only; does not alter preregistered gates.
        if name == 'n9_cfl005':
            x = history['reference_position_m'][:, 0]
            volume = history['volume_m3']
            sigma = history['stress_Pa'][-1]
            ref = history['reference_final_stress_Pa']
            error2 = volume*np.sum((sigma-ref)**2, axis=(1, 2))
            axial_fraction = float(np.sum(volume*(sigma[:, 0, 0]-ref[:, 0, 0])**2)/error2.sum())
            for coordinate in np.unique(x):
                on = x == coordinate
                profiles.append(dict(x_m=float(coordinate),
                    fraction_of_total_squared_stress_error=float(error2[on].sum()/error2.sum()),
                    actual_weighted_plane_sigma_xx_Pa=float(np.average(sigma[on, 0, 0], weights=volume[on])),
                    reference_weighted_plane_sigma_xx_Pa=float(np.average(ref[on, 0, 0], weights=volume[on]))))
    try:
        restore(case/'checkpoint.npz')
    except ValueError:
        pass
    else:
        raise AssertionError('incomplete benchmark restart was incorrectly authorized')
    last = report['records'][-1]
    residual = (np.array(last['linear_momentum_kg_m_s'])
                -np.array(last['constraint_impulse_N_s'])-np.array(last['body_impulse_N_s']))
    relative = float(np.linalg.norm(residual)/np.linalg.norm(last['body_impulse_N_s']))
    assert relative < 1e-10, (name, relative)
    errors.append(relative)
print(json.dumps(dict(status='PASS', scope='artifact integrity and impulse accounting only',
                     cases=len(names), max_momentum_relative_residual=max(errors),
                     production_authorized=False, scientific_status='UNMODIFIED_FAIL_ENGINEERING_TARGETS',
                     posthoc_localization_not_new_acceptance=dict(n9_axial_squared_error_fraction=axial_fraction,
                     n9_plane_profile=profiles)), indent=2))
