"""Inherited geometric guards and the saved event-location high-probe chord."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import afsd_event_blocked_box_20261010 as proof
import test_event_blocked_path_20261010 as legacy


class BoxOnlyTests(legacy.BlockedPrefixTests):
    """All thirteen independent legacy geometric cases use the new audit."""
    def audit(self,end=.49,**replace):
        args=dict(x0=self.x0,x1=self.x1,i=self.i,j=self.j,pose0=self.pose,pose1=self.pose,
            p=self.p,particle_ids=self.ids,end_fraction=end)
        args.update(replace)
        return proof.audit_blocked_prefix(**args)


class SavedBoxTests(legacy.SavedBlockedPrefixTests):
    """All four SHA-bound historical failed-event cases use the new audit."""
    def audit(self,end,**kw):
        return proof.audit_blocked_prefix(self.x0,self.x1,np.array([1333]),np.array([1370]),
            self.pose0,self.pose1,self.p,self.ids,end,**kw)


class ProofPreservationTests(unittest.TestCase):
    def test_original_math_domain_and_arithmetic_helpers_are_identical(self):
        old_path=ROOT/'scripts/afsd_event_blocked_path_20261010.py'
        self.assertEqual(hashlib.sha256(old_path.read_bytes()).hexdigest(),proof.FROZEN_PARENT_SHA256)
        old_tree=ast.parse(old_path.read_text())
        new_tree=ast.parse((ROOT/'scripts/afsd_event_blocked_box_20261010.py').read_text())
        old=next(n for n in old_tree.body if isinstance(n,ast.FunctionDef))
        new=next(n for n in new_tree.body if isinstance(n,ast.FunctionDef))
        # Domain checks, strict sampling, fixed-s witness, roundoff and
        # curvature are exactly the frozen implementations.
        for name in ('stop','local_at','query','strict_at','box_witnesses'):
            before=next(n for n in old.body if isinstance(n,ast.FunctionDef) and n.name==name)
            after=next(n for n in new.body if isinstance(n,ast.FunctionDef) and n.name==name)
            self.assertEqual(ast.dump(before),ast.dump(after),name)
        # Input/budget validation, rotation, speed and margin construction.
        def before_report(fn):
            nodes=[]
            for node in fn.body[1:]:
                if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='report' for t in node.targets):
                    break
                nodes.append(ast.dump(node))
            return nodes
        self.assertEqual(before_report(old),before_report(new))
        old_try=next(n for n in old.body if isinstance(n,ast.Try))
        new_try=next(n for n in new.body if isinstance(n,ast.Try))
        old_while=next(n for n in old_try.body if isinstance(n,ast.While))
        new_while=next(n for n in new_try.body if isinstance(n,ast.While))
        # Full-prefix minimum-length/degeneration test and strict endpoints.
        self.assertEqual([ast.dump(n) for n in old_try.body[:old_try.body.index(old_while)]],
            [ast.dump(n) for n in new_try.body[:new_try.body.index(new_while)]])
        # The box-only loop retains the entire old box proof, then the old
        # depth-limited binary subdivision, without an additional gate.
        delta=next(k for k,n in enumerate(old_while.body) if isinstance(n,ast.Assign)
            and any(isinstance(t,ast.Name) and t.id=='delta' for t in n.targets))
        self.assertEqual([ast.dump(n) for n in old_while.body[:delta]],
            [ast.dump(n) for n in new_while.body[:delta]])
        self.assertEqual(ast.dump(old_while.body[-1]),ast.dump(new_while.body[-1]))
        self.assertEqual([ast.dump(n) for n in old_try.body[old_try.body.index(old_while)+1:]],
            [ast.dump(n) for n in new_try.body[new_try.body.index(new_while)+1:]])
        self.assertEqual([ast.dump(n) for n in old_try.handlers],[ast.dump(n) for n in new_try.handlers])


class SavedHighProbeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        stage=ROOT/'outputs/gpu_acceleration_20261007/event_secant_sensitivity_20261010_v1/sensitivity_one/coarse/event_probe_08'
        expected={
            'visible_graph.npz':'0328b03d5bcfc739b6951fa08093b34f90d0e7ce68911435267c9cbdfc2abeef',
            'motion.json':'f3e4eb59e323a2c891cde84ec33de4603aad8be024133cf17b9a5bd94383a0a6',
            'incoming_blocked_prefix.json':'dbd47d05dd969d84eb3d9a32f8e28baa93200383c6877c964d13e660e9383a84',
        }
        for name,digest in expected.items():
            if hashlib.sha256((stage/name).read_bytes()).hexdigest()!=digest:
                raise AssertionError('Saved high-probe evidence changed: '+name)
        with np.load(stage/'visible_graph.npz',allow_pickle=False) as z:
            cls.x0,cls.x1,cls.ids=z['position0'],z['position1'],z['particle_id']
        cls.i=np.flatnonzero(cls.ids==1334);cls.j=np.flatnonzero(cls.ids==1371)
        motion=json.loads((stage/'motion.json').read_text())
        cls.pose0=dict(origin=motion['origin0_m'],rotation=motion['rotation0'])
        cls.pose1=dict(origin=motion['origin1_m'],rotation=motion['rotation1'])
        cls.old_report=json.loads((stage/'incoming_blocked_prefix.json').read_text())
        cls.end=cls.old_report['prefix_end_fraction']
        cfg=json.loads((ROOT/'configs/deposition_fill_then_track_20261007.json').read_text())
        cls.p=json.loads((ROOT/cfg['base_parameters']).read_text())['parameters']
        cls.p.update({k:v['value'] for k,v in cfg['parameter_overrides'].items()})

    def audit(self,end=None,**kw):
        return proof.audit_blocked_prefix(self.x0,self.x1,self.i,self.j,
            self.pose0,self.pose1,self.p,self.ids,self.end if end is None else end,**kw)

    def test_same_failed_high_chord_is_certified_with_original_box_and_budget(self):
        self.assertEqual(self.old_report['status'],proof.UNRESOLVED)
        self.assertEqual(self.old_report['failure_kind'],'MAXIMUM_ANALYTIC_BATCH_QUERIES')
        self.assertEqual(self.old_report['analytic_batch_queries'],256)
        report=self.audit()
        self.assertEqual(report['input_sha256'],self.old_report['input_sha256'])
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertEqual(report['failure_kind'],None)
        self.assertLessEqual(report['analytic_batch_queries'],256)
        self.assertLessEqual(report['deepest_subdivision'],32)
        self.assertEqual(report['contracted_solid_queries'],0)
        self.assertEqual(report['analytic_batch_queries'],
            report['strict_geometry_batch_queries']+report['fixed_witness_box_queries'])
        self.assertEqual(report['partition_queries'],report['strict_geometry_batch_queries'])
        self.assertEqual({r['kind'] for r in report['query_log']},{'strict','fixed_witness_box'})
        self.assertEqual({r['method'] for r in report['certificates']},
            {'FIXED_SEGMENT_PARAMETER_CURVATURE_BOX'})
        intervals=sorted((r['start_fraction'],r['end_fraction']) for r in report['certificates'])
        self.assertEqual(intervals[0][0],0.)
        self.assertEqual(intervals[-1][1],self.end)
        self.assertTrue(all(a[1]==b[0] for a,b in zip(intervals,intervals[1:])))
        self.assertEqual(report['physical_solver_evaluations'],0)
        self.assertEqual(report['edge_deletions'],0)
        json.dumps(report)

    def test_same_high_chord_through_visible_endpoint_is_rejected(self):
        report=self.audit(1.)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['failure_kind'],'INCOMING_NOT_STRICTLY_BLOCKED_AT_SAMPLE')
        self.assertTrue(report['visible_samples'])
        self.assertEqual(report['visible_samples'][0]['fraction'],1.)
        self.assertEqual(report['analytic_batch_queries'],2)

    def test_high_chord_budget_or_depth_exhaustion_never_becomes_pass(self):
        for kw,reason in ((dict(maximum_queries=223),'MAXIMUM_ANALYTIC_BATCH_QUERIES'),
                (dict(maximum_depth=17),'MAXIMUM_SUBDIVISION_DEPTH')):
            with self.subTest(kw=kw):
                report=self.audit(**kw)
                self.assertEqual(report['status'],proof.UNRESOLVED)
                self.assertEqual(report['failure_kind'],reason)
                self.assertTrue(report['unresolved_intervals'])
                self.assertEqual(report['contracted_solid_queries'],0)
                self.assertLessEqual(report['analytic_batch_queries'],kw.get('maximum_queries',256))


if __name__=='__main__':
    unittest.main()
