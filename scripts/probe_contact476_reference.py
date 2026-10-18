"""Independent frozen two-face cone reference; NOT_REPAIR and no AFSD advancement.

The objective is b.f + f.H.f/2 on two force disks, H=c*A.T*M*A.
For each two-dimensional block, solve (H_aa + lambda*I) f_a=-b_eff;
lambda=0 for an interior solution, otherwise bisect the disk-radius equation.
This is a direct convex KKT derivation, not a call into contact_cone/contact_kick.
The 80-digit certificate concerns the stored decimal input, not unknown true data.
"""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import time


def reference(context_path):
    start = time.perf_counter()
    raw = Path(context_path).read_bytes()
    context = json.loads(raw.decode('utf-8'), parse_float=Decimal)
    local = context['failed_local_context']
    cone = local['failed_cone']
    with localcontext() as arithmetic:
        arithmetic.prec = 80
        zero = Decimal(0)
        coefficient = cone['coefficient']
        bases = cone['tangent_bases']
        mobility = cone['mobility']
        capacities = cone['limits_N']
        b = sum(cone['base_tangent'], [])
        candidate = sum(cone['force_tangent_N'], [])
        if len(capacities) != 2 or len(b) != 4 or len(candidate) != 4:
            raise ValueError('this independent reference is limited to two faces')
        if len(bases) != 2 or any(len(a) != 3 or any(len(row) != 2 for row in a) for a in bases):
            raise ValueError('two 3x2 tangent bases required')
        if len(mobility) != 3 or any(len(row) != 3 for row in mobility):
            raise ValueError('3x3 mobility required')
        # The recorded candidate has diagonal positive mobility. Refuse a new
        # class rather than silently claiming a general semidefinite solver.
        if any(mobility[i][j] != 0 for i in range(3) for j in range(3) if i != j):
            raise ValueError('reference scope: positive diagonal mobility only')
        if coefficient <= 0 or any(r <= 0 for r in capacities) or any(mobility[i][i] <= 0 for i in range(3)):
            raise ValueError('positive coefficient, capacities and mobility required')
        scalar_input = [coefficient, *capacities, *b, *candidate]
        scalar_input += [v for matrix in bases for row in matrix for v in row]
        scalar_input += [v for row in mobility for v in row]
        if any(not Decimal(v).is_finite() for v in scalar_input):
            raise ValueError('nonfinite frozen input')

        norm = lambda values: sum((v*v for v in values), zero).sqrt()
        a_matrix = [[bases[j//2][i][j % 2] for j in range(4)] for i in range(3)]
        hessian = [[coefficient*sum(
            (a_matrix[k][i]*mobility[k][ell]*a_matrix[ell][j]
             for k in range(3) for ell in range(3)), zero)
            for j in range(4)] for i in range(4)]
        force = [zero]*4
        multipliers = [zero]*2
        sweep_change = None
        for sweep in range(12):
            previous = force[:]
            for face in range(2):
                indices = [2*face, 2*face+1]
                effective_b = [b[i]+sum((hessian[i][j]*force[j]
                    for j in range(4) if j not in indices), zero) for i in indices]
                h00, h01 = (hessian[indices[0]][j] for j in indices)
                h10, h11 = (hessian[indices[1]][j] for j in indices)

                def block_force(multiplier):
                    determinant = (h00+multiplier)*(h11+multiplier)-h01*h10
                    if determinant <= 0:
                        raise ValueError('reference requires positive definite diagonal blocks')
                    return [(-(h11+multiplier)*effective_b[0]+h01*effective_b[1])/determinant,
                            (h10*effective_b[0]-(h00+multiplier)*effective_b[1])/determinant]

                unconstrained = block_force(zero)
                if norm(unconstrained) <= capacities[face]:
                    updated = unconstrained
                    multiplier = zero
                else:
                    lower = zero
                    upper = norm(effective_b)/capacities[face]
                    if norm(block_force(upper)) > capacities[face]:
                        raise ArithmeticError('independent multiplier bracket failed')
                    for _ in range(220):
                        middle = (lower+upper)/2
                        if norm(block_force(middle)) > capacities[face]:
                            lower = middle
                        else:
                            upper = middle
                    multiplier = upper
                    updated = block_force(multiplier)
                for index, value in zip(indices, updated):
                    force[index] = value
                multipliers[face] = multiplier
            sweep_change = norm([force[i]-previous[i] for i in range(4)])
            if sweep_change < Decimal('1e-60'):
                break

        def gradient(values):
            return [b[i]+sum((hessian[i][j]*values[j] for j in range(4)), zero) for i in range(4)]

        def objective(values):
            return sum((b[i]*values[i] for i in range(4)), zero)+sum(
                (values[i]*hessian[i][j]*values[j]/2 for i in range(4) for j in range(4)), zero)

        def gap(values):
            g = gradient(values)
            return sum((g[i]*values[i] for i in range(4)), zero)+sum(
                (capacities[k]*norm(g[2*k:2*k+2]) for k in range(2)), zero)

        force_scale = max(capacities)
        mobility_norm = max(mobility[i][i] for i in range(3))
        velocity_scale = max(norm(b), coefficient*mobility_norm*force_scale)
        g = gradient(force)
        stationarity = norm([g[i]+multipliers[i//2]*force[i] for i in range(4)])
        reference_gap = gap(force)/(force_scale*velocity_scale)
        reference_disk_excess = max(norm(force[2*k:2*k+2])-capacities[k] for k in range(2))
        complementarity = max(abs(multipliers[k]*(norm(force[2*k:2*k+2])-capacities[k])) for k in range(2))
        # These are high-precision numerical certificate thresholds, not proposed
        # replacements for any frozen production/contact acceptance tolerance.
        checks = {
            'iterate_change_N': sweep_change < Decimal('1e-60'),
            'disk_feasibility_N': reference_disk_excess <= Decimal('1e-60'),
            'stationarity_base_units': stationarity <= Decimal('1e-60'),
            'complementarity_objective_units': complementarity <= Decimal('1e-60'),
            'signed_scaled_gap': abs(reference_gap) <= Decimal('1e-60'),
            'nonnegative_multipliers': all(value >= 0 for value in multipliers),
        }
        resultant_error = [sum((a_matrix[i][j]*(candidate[j]-force[j])
            for j in range(4)), zero) for i in range(3)]
        # In midpoint(), coefficient=h/2. This is only the frozen local
        # subinterval; it is NOT a bound on the complete kick/AFSD trajectory.
        interval = 2*coefficient
        delta_velocity = [interval*sum((mobility[i][j]*resultant_error[j]
            for j in range(3)), zero) for i in range(3)]
        report = {
            'classification': 'FROZEN_LOCAL_DIAGNOSTIC_NOT_REPAIR',
            'status': 'REFERENCE_CERTIFIED' if all(checks.values()) else 'REFERENCE_FAILED',
            'real_case_status': 'FAILED_UNCHANGED',
            'production_authorized': False,
            'method': 'independent convex two-disk block minimization with scalar KKT multiplier bisection',
            'reference_source': 'Direct derivation of KKT conditions for stored convex quadratic; no tested solver calls.',
            'precision_decimal_digits': 80,
            'maximum_sweeps': 12,
            'multiplier_bisections_per_active_block': 220,
            'certificate_threshold': '1e-60 in the explicitly named units/scales; diagnostic only',
            'sweeps': sweep+1,
            'checks': checks,
            'particle_id': int(local['particle_id']),
            'features': [face['feature'] for face in local['faces']],
            'reference_force_tangent_N': [[str(v) for v in force[:2]], [str(v) for v in force[2:]]],
            'candidate_force_tangent_N': [[str(v) for v in candidate[:2]], [str(v) for v in candidate[2:]]],
            'reference_kkt_stationarity_absolute_base_units': str(stationarity),
            'reference_scaled_gap': str(reference_gap),
            'reference_max_disk_excess_N': str(reference_disk_excess),
            'reference_complementarity_absolute_objective_units': str(complementarity),
            'candidate_scaled_gap_high_precision': str(gap(candidate)/(force_scale*velocity_scale)),
            'candidate_objective_suboptimality_W': str(objective(candidate)-objective(force)),
            'candidate_objective_suboptimality_scaled': str((objective(candidate)-objective(force))/(force_scale*velocity_scale)),
            'candidate_max_disk_excess_N': str(max(norm(candidate[2*k:2*k+2])-capacities[k] for k in range(2))),
            'force_vector_error_N': str(norm([candidate[i]-force[i] for i in range(4)])),
            'resultant_force_error_N': str(norm(resultant_error)),
            'local_midpoint_interval_h_s': str(interval),
            'conditional_frozen_local_endpoint_velocity_difference_m_s': str(norm(delta_velocity)),
            'conditional_frozen_local_endpoint_velocity_difference_upper_bound_m_s': str(interval*mobility_norm*norm(resultant_error)),
            'candidate_reported_scaled_residual': str(cone['scaled_residual']),
            'candidate_reported_internal_target': str(cone['target']),
            'context_file': str(Path(context_path).resolve()),
            'context_sha256': hashlib.sha256(raw).hexdigest(),
            'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'wall_seconds': time.perf_counter()-start,
            'interpretation': [
                'The recorded global certificate misses its target; this does not establish a macroscopic physical failure.',
                'A near-zero-slip face can have a large force error while its tiny local impulse has a roundoff-scale effect.',
                'The local interval and velocity bound do not certify the full contact event or justify accepting the rejected state.',
                'No tolerance, force, heat, trajectory, checkpoint, source law or true-run status is repaired here.',
                'Objective units W apply because this recorded coefficient is the midpoint time coefficient and base is velocity.',
            ],
        }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--context', required=True)
    parser.add_argument('--output', required=True, help='new JSON path; existing file is never overwritten')
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(output)
    report = reference(args.context)
    output.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+'\n'
    with output.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(text)
    print(text, end='')
    return 0 if report['status'] == 'REFERENCE_CERTIFIED' else 1


if __name__ == '__main__':
    raise SystemExit(main())
