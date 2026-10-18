"""Independent algebraic counterexamples for route selection, not a GPU run."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FE = ROOT / 'third_party/route_audit_20261004/FE-SPH-GPU'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    dest = args.output.resolve()
    if not dest.is_relative_to((ROOT/'outputs').resolve()) or dest.exists():
        raise ValueError('new file inside outputs required')
    sources = [FE/'src/interactions_gpu.cu', FE/'src/types.h',
               FE/'glm/detail/type_vec3.inl', ROOT/'src/afsd_backend/pysph/material.py',
               ROOT/'outputs/contact_material_search_20261004/root/resolved_config.json']
    before = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sources}
    if 'vr = v - v * normal;' not in sources[0].read_text():
        raise ValueError('reviewed source expression is absent')
    rows = []
    for n in (np.array([1/math.sqrt(2), 1/math.sqrt(2), 0.]), np.array([0., 0., -1.])):
        v = n.copy()  # Purely normal motion has exactly zero tangential component.
        candidate = v-v*n  # Reviewed GLM vector product is componentwise.
        orthogonal = v-np.dot(v, n)*n
        if np.linalg.norm(orthogonal) > 1e-14 or abs(float(np.dot(candidate, n))) < .1:
            raise ArithmeticError('independent projection counterexample inconsistent')
        rows.append(dict(unit_normal=n.tolist(), relative_velocity=v.tolist(),
                         reviewed_expression_result=candidate.tolist(),
                         correct_tangential_projection=orthogonal.tolist(),
                         spurious_normal_component=float(np.dot(candidate, n))))
    p = json.loads(sources[-1].read_text())['parameters']
    ts = (p['temperature']-p['room'])/(p['melt']-p['room'])
    coefficient = p['b1']*ts**p['b2']
    threshold = p['reference_rate']*math.exp(-1/coefficient)
    prefactor_rows = [dict(rate_s_inv=r, unregularized_rate_factor=1+coefficient*math.log(r/p['reference_rate']))
                     for r in (.1, .5, 1.)]
    after = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in before}
    if before != after:
        raise RuntimeError('protected input changed')
    result = dict(classification='SOURCE_EXPRESSION_COUNTEREXAMPLES_NOT_UPSTREAM_EXECUTION',
                  contact_projection=dict(status='COUNTEREXAMPLE_CONFIRMED', examples=rows,
                                          source='FE-SPH-GPU/src/interactions_gpu.cu:933,1001',
                                          scope='Independent algebra only; no CUDA build or trajectory prediction.'),
                  FKS_low_rate=dict(temperature_K=p['temperature'],
                                    zero_of_unregularized_rate_factor_s_inv=threshold,
                                    current_declared_floor_s_inv=p['rate_floor'],
                                    samples=prefactor_rows,
                                    scope='Removing the existing floor is not a valid slow-flow model; production law unchanged.'),
                  input_sha256=before, inputs_unchanged=True,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  new_physical_time_s=0., AFSD_validation=False)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(dest), projection='COUNTEREXAMPLE_CONFIRMED',
                          unregularized_rate_factor_zero_s_inv=threshold)))


if __name__ == '__main__':
    main()
