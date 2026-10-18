"""Retarget original fixtures plus independent supplemental-bound regressions."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import types
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import visible_support_anchored_20261018 as proof
import visible_support_prism_20261010 as previous


# Retarget every unchanged old fixture in an isolated test module. Production
# globals are not rebound, and old files are never written.
fixture_path=ROOT/'tests/review/test_visible_support_prism_20261010.py'
fixture_source=fixture_path.read_text().replace(
    'import visible_support_prism_20261010 as proof',
    'import visible_support_anchored_20261018 as proof')
fixtures=types.ModuleType('anchored_inherited_prism_fixtures')
fixtures.__file__=str(fixture_path)
exec(compile(fixture_source,str(fixture_path),'exec'),fixtures.__dict__)


class SupplementalBoundTests(unittest.TestCase):
    def test_original_geometry_functions_and_loop_are_preserved(self):
        def functions(path):
            tree=ast.parse(path.read_text())
            audit=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='audit_visible_support_path')
            return {n.name:ast.dump(n) for n in audit.body if isinstance(n,ast.FunctionDef)},audit
        old,oa=functions(ROOT/'scripts/visible_support_prism_20261010.py')
        new,na=functions(ROOT/'scripts/visible_support_anchored_20261018.py')
        for name in ('query','local_at','strict_at','unresolved'):
            self.assertEqual(old[name],new[name],name)
        self.assertEqual(ast.dump(next(n for n in ast.walk(oa) if isinstance(n,ast.While))),
                         ast.dump(next(n for n in ast.walk(na) if isinstance(n,ast.While))))
        self.assertEqual(hashlib.sha256((ROOT/'scripts/visible_support_prism_20261010.py').read_bytes()).hexdigest(),
            'a71f9eb392d06674a0cae748b7238cbedaa187684bcb66566cc96d21b1833bab')

    def test_matrix_envelope_matches_analytic_inequalities(self):
        for angle in (1e-6,1.85e-4,.2,1.,np.pi/2):
            R0=np.eye(3);c,s=np.cos(angle),np.sin(angle)
            R1=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
            first=np.array([[1.,2.,3.],[4.,3.,2.]])+1e5
            last=first+[.1,-.2,.3];origin=np.full(3,1e5)
            selected,report=proof._anchored_position_margin(first,last,origin,origin,R0,R1,
                angle,1e5,1.,dict(per_particle_curvature_upper_m=[1.,1.],original_displacement_roundoff_margin_m=1e-14))
            self.assertTrue(report['available'])
            self.assertEqual(report['roundoff_factor'],8192*np.finfo(float).eps)
            self.assertGreater(report['trigonometric_error_m'],0.)
            self.assertGreater(selected,2*report['uniform_position_evaluation_error_m'])
            for f in np.linspace(0.,1.,33):
                c,s=np.cos(angle*f),np.sin(angle*f)
                matrix=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
                self.assertLessEqual(np.linalg.norm(matrix,ord='fro'),report['matrix_norm_upper'])

    def test_large_world_translation_remains_in_error_scale(self):
        def bound(shift):
            first=np.array([[.1,.2,.3],[.4,.3,.2]])+shift;last=first+[1e-4,0.,0.]
            origin=np.full(3,shift);angle=.01;c,s=np.cos(angle),np.sin(angle)
            return proof._anchored_position_margin(first,last,origin,origin,np.eye(3),
                np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]),angle,max(1.,abs(shift)),1.,
                dict(per_particle_curvature_upper_m=[1.,1.],original_displacement_roundoff_margin_m=1e-14))[1]
        self.assertGreater(bound(1e6)['arithmetic_error_m'],10000*bound(0.)['arithmetic_error_m'])

    def test_outside_inherited_domain_has_no_supplement(self):
        x=np.array([[.1,.2,.3],[.4,.3,.2]])
        for angle in (0.,5e-7,2.):
            margin,r=proof._anchored_position_margin(x,x,np.zeros(3),np.zeros(3),np.eye(3),np.eye(3),
                angle,1.,.01,dict(per_particle_curvature_upper_m=[1.,1.]))
            self.assertEqual(margin,.01);self.assertFalse(r['available'])

    def test_saved_after_event_exact_chord_is_certified(self):
        base=ROOT/'outputs/short_bar_event_continue_20261015_v4/run/piece_000002_attempt_000001/fine_1/after_event/proposal'
        hashes={name:hashlib.sha256((base/name).read_bytes()).hexdigest() for name in
            ('visible_graph.npz','motion.json','visible_path_audit.json','NOT_RESTARTABLE_candidate.npz')}
        old=json.loads((base/'visible_path_audit.json').read_text())
        self.assertEqual(old['status'],proof.UNRESOLVED)
        cfg=json.loads((ROOT/'configs/deposition_fill_then_track_20261007.json').read_text())
        p=json.loads((ROOT/cfg['base_parameters']).read_text())['parameters']
        p.update({k:v['value'] for k,v in cfg['parameter_overrides'].items()})
        motion=json.loads((base/'motion.json').read_text())
        with np.load(base/'visible_graph.npz',allow_pickle=False) as z:
            report=proof.audit_visible_support_path(z['position0'],z['position1'],z['i'],z['j'],
                dict(origin=motion['origin0_m'],rotation=motion['rotation0']),
                dict(origin=motion['origin1_m'],rotation=motion['rotation1']),p,z['particle_id'])
        self.assertEqual(report['input_sha256'],old['input_sha256'])
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertLessEqual(report['analytic_batch_queries'],256)
        self.assertFalse(report['found_edges'] or report['degenerate_edges'] or report['unresolved_intervals'])
        self.assertTrue(any(w['endpoint_enclosure_method']=='ANGLE_AWARE_ENDPOINT_ANCHORED_CURVATURE_PRISM'
            for row in report['certificates'] for w in row.get('witnesses',[])))
        for name,value in hashes.items():self.assertEqual(hashlib.sha256((base/name).read_bytes()).hexdigest(),value)


def load_tests(loader,tests,pattern):
    suite=loader.loadTestsFromModule(fixtures)
    suite.addTests(loader.loadTestsFromTestCase(SupplementalBoundTests))
    return suite


if __name__=='__main__':unittest.main()
