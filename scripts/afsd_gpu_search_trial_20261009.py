"""Versioned GPU trial with an explicit auxiliary material-tangent provider.

The residual/evaluate, thermal solve, physical material update, state advance,
contact, geometry and prescribed machine motion are inherited unchanged.
Only tangent imports the new finite-symmetric-stress derivative provider.
No historical module global is rebound. This module alone authorizes no run.
"""
import ast
import copy
import inspect
import textwrap
import numpy as np
import afsd_gpu_trial_20261007 as parent
from afsd_search_tensor_20261009 import static_adaptation_audit, SEARCH_PROTOCOL

ROOT, CLASSIFICATION = parent.ROOT, parent.CLASSIFICATION
DEFAULT_CONFIG, NUMERICAL_PROTOCOL = parent.DEFAULT_CONFIG, parent.NUMERICAL_PROTOCOL
COMPUTE_PROTOCOL = parent.COMPUTE_PROTOCOL
require, identity_digest = parent.require, parent.identity_digest
validate_state, independent_initial_state = parent.validate_state, parent.independent_initial_state
load_configuration, advance_candidate = parent.load_configuration, parent.advance_candidate


class ProcessTrial(parent.ProcessTrial):
    def tangent(self, result):
        """Neighbor-local assembly + local contact secant; NOT full Jacobian."""
        from afsd_backend.pysph.quasistatic_tangent import (gradient_map,
                                          assemble_frozen_force_tangent)
        from afsd_search_tensor_20261009 import material_tangent
        self._require_interval(result['dt_s'])
        graph=result['visibility_graph']
        for key in ('i','j','g'):
            if not np.array_equal(result[key],graph[key]):
                raise ValueError('tangent must use the exact residual visible graph')
        p, r = self.p, result
        dt = r['dt_s']
        D = gradient_map(r['volume'], r['i'], r['j'], r['g'], r['correction'])
        material = material_tangent(self.s0, self.pa.ep, r['gradient'],
            r['material_temperature'], dt, p['shear'], self.law, material=r['material'])
        J = assemble_frozen_force_tangent(D, r['volume'], material['tangent_Pa_s'],
            self.pa.rho, p['rho0'], self.bulk, dt, r['displacement'])
        blocks = np.zeros((len(self.x0), 3, 3))
        # Batch perturb independent particle contacts, six contact calls total.
        eps = 1e-9
        for axis in range(3):
            plus = r['position'].copy(); minus = plus.copy()
            plus[:, axis] += eps; minus[:, axis] -= eps
            args = (self.motion.pose(self.state['time_s']), self.motion.pose(self.state['time_s']+dt),
                    dt, self.state['contact_history'])
            a = self.contact.evaluate(self.x0, plus, *args, check_path=False)
            b = self.contact.evaluate(self.x0, minus, *args, check_path=False)
            blocks[:, :, axis] = (a['force_N']-b['force_N'])/(2*eps)
        for particle, block in enumerate(blocks):
            rows = slice(3*particle, 3*particle+3)
            J[rows, rows] += block
        return J


def source_hashes(config_path=DEFAULT_CONFIG):
    mapping = parent.source_hashes(config_path)
    import hashlib
    for name in ('afsd_search_tensor_20261009.py', 'afsd_gpu_search_trial_20261009.py'):
        path = ROOT/'scripts'/name
        mapping[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return mapping


def operator_reuse_audit(actual_tree=None):
    inherited = parent.operator_reuse_audit()
    reference = ast.parse(textwrap.dedent(inspect.getsource(parent.ProcessTrial.tangent))).body[0]
    original_import = reference.body[1]
    require(isinstance(original_import, ast.ImportFrom)
            and original_import.module == 'afsd_backend.pysph.quasistatic_tangent',
            'unexpected original tangent import')
    original_import.names = [n for n in original_import.names if n.name != 'material_tangent']
    reference.body.insert(2, ast.ImportFrom(module='afsd_search_tensor_20261009',
                          names=[ast.alias(name='material_tangent', asname=None)], level=0))
    actual = (ast.parse(textwrap.dedent(inspect.getsource(ProcessTrial.tangent))).body[0]
              if actual_tree is None else copy.deepcopy(actual_tree))
    require(ast.dump(reference, include_attributes=False) == ast.dump(actual, include_attributes=False),
            'tangent changed beyond explicit material search provider')
    for name in ('evaluate', 'prescribed', '__init__', '_require_interval'):
        require(getattr(ProcessTrial, name) is getattr(parent.ProcessTrial, name),
                'inherited physical trial provider changed: '+name)
    require(advance_candidate is parent.advance_candidate and validate_state is parent.validate_state,
            'state advance or validation changed')
    derivative = static_adaptation_audit()
    return dict(status='PASS_EXPLICIT_SEARCH_TANGENT_PROVIDER_ADAPTATION',
                inherited_audit=inherited, derivative_audit=derivative,
                endpoint_evaluate_inherited_unchanged=True, physical_gates_changed=False,
                runtime_original_global_rebinding=False, stress_projected=False,
                continuation_authorized=False, cuda_executed=False)
