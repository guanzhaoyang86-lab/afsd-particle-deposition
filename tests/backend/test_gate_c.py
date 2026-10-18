import tempfile
import unittest
from pathlib import Path
import numpy as np
from afsd_backend.pysph.config_adapter import load_settings
from afsd_backend.pysph.particles import create_geometry, identity_map, append_feedstock
from afsd_backend.pysph.output import checkpoint, restore, vtk
from afsd_backend.pysph.geometry import as_rs
from afsd_scaffold.vtk_legacy import read_point_cloud
from afsd_scaffold.config import ConfigurationError

ROOT = Path(__file__).resolve().parents[2]


class GateC(unittest.TestCase):
    def setUp(self):
        self.config, self.p, _ = load_settings(ROOT/'configs/debug/starved.json')

    def test_three_cases_and_no_paper_fallback(self):
        for name, feed in [('starved', .00106), ('copper-core', .00212), ('overfed', .00424)]:
            _, p, _ = load_settings(ROOT/f'configs/debug/{name}.json')
            self.assertEqual(p['feed'], feed)
            self.assertAlmostEqual(p['omega'], 10*np.pi)
            self.assertEqual(p['layer'], .0015)
        with self.assertRaises(ConfigurationError): load_settings(ROOT/'configs/cases/starved.json')

    def test_id_reorder_restart_inflow_and_vtk(self):
        pa = create_geometry(self.p); before = identity_map(pa)
        pa.S[:] = np.arange(pa.S.size, dtype=float)
        stress_by_id = dict(zip(pa.particle_id.tolist(), pa.S.reshape(-1, 9).tolist()))
        order = np.random.default_rng(73).permutation(pa.get_number_of_particles()).astype(np.uint32)
        pa = pa.extract_particles(order)
        self.assertEqual(identity_map(pa), before)
        self.assertEqual(dict(zip(pa.particle_id.tolist(), pa.S.reshape(-1, 9).tolist())), stress_by_id)
        with tempfile.TemporaryDirectory() as d:
            checkpoint(Path(d)/'state.npz', pa, {'t': .1, 'next_id': 10000})
            pb, state = restore(Path(d)/'state.npz')
            self.assertEqual(identity_map(pb), before); self.assertEqual(state['next_id'], 10000)
            vtk(Path(d)/'state.vtk', pb, self.p, 68e9, .1)
            out = read_point_cloud(Path(d)/'state.vtk')
            self.assertEqual(out['scalars']['particle_id'], pb.particle_id.tolist())
            self.assertEqual(out['vectors']['initial_position_m'], list(map(tuple, np.column_stack((pb.x0, pb.y0, pb.z0)).tolist())))
        append_feedstock(pa, [[0., 0., .005]], self.p, 10000)
        after = identity_map(pa)
        self.assertEqual({k: after[k] for k in before}, before)
        self.assertIn(10000, after)
        with self.assertRaises(ValueError): append_feedstock(pa, [[0., 0., .006]], self.p, 10000)

    def test_local_tags_and_frames(self):
        pa = create_geometry(self.p)
        self.assertTrue(np.any(pa.tracer_type == 1)); self.assertTrue(np.any(pa.tracer_type == 2))
        self.assertTrue(np.all((pa.z0[pa.tracer_type > 0] >= .003) & (pa.z0[pa.tracer_type > 0] <= .004)))
        self.assertGreater(np.count_nonzero(pa.wholebar_tag), np.count_nonzero(pa.tracer_type))
        np.testing.assert_equal(pa.material_id, 1)
        np.testing.assert_equal(as_rs([-1, 0, 1], 1, 1), [1, 0, -1])
        np.testing.assert_equal(as_rs([-1, 0, 1], -1, 1), [-1, 0, 1])

    def test_unevaluated_initial_rates_are_absent(self):
        pa = create_geometry(self.p)
        with tempfile.TemporaryDirectory() as d:
            fields = vtk(Path(d)/'initial.vtk', pa, self.p, 68e9, 0.)
            self.assertFalse(any('rate' in field for field in fields))
            state = read_point_cloud(Path(d)/'initial.vtk')
            self.assertFalse(any('rate' in field for field in state['scalars']))

    def test_id_overflow_fails_instead_of_wrapping(self):
        from afsd_backend.pysph.particles import make_array
        with self.assertRaises(OverflowError):
            make_array([[0, 0, .005]], np.ones(1, int), self.p, 2**63-1)
