"""Read saved contact/force evidence only; no physical step or state edits."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import run_gpu_nonlinear_resume_20261011 as resume


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    diagnostic = read(Path(__file__).with_name('diagnostic.json'))
    rows = []
    for source, frame in zip(diagnostic['sources'], diagnostic['records']):
        checkpoint = ROOT / source['checkpoint']
        assert sha(checkpoint) == source['sha256']
        pa, state = resume.load(checkpoint)
        proposal = checkpoint.parent / 'proposal'
        jp = proposal / 'SOLVED_UNCOMMITTED_candidate.json'
        npz = proposal / 'SOLVED_UNCOMMITTED_candidate.npz'
        motion = read(proposal / 'motion.json')
        evidence = read(jp)
        contact = evidence['contact']
        energy = evidence['energy']
        with np.load(npz, allow_pickle=False) as saved:
            ids = saved['particle_id'].copy()
            pos = saved['position_m'].copy()
            total = saved['force_N'].copy()
            internal = saved['internal_force_N'].copy()
            contact_force = saved['contact_force_N'].copy()
            velocity = saved['velocity_m_s'].copy()
            assert np.array_equal(ids, pa.particle_id)
            assert np.array_equal(pos, np.column_stack((pa.x, pa.y, pa.z)))
            assert np.array_equal(saved['stress_Pa'], pa.S.reshape(-1, 3, 3))
        assert abs(motion['t1_s'] - state['time_s']) < 1e-15
        free = np.asarray(state['free'], bool)
        actuator = np.asarray(pa.actuator, bool)
        reaction = np.asarray(energy['prescribed_reaction_N'])
        normal = np.asarray(contact['normal_force_by_feature_N'])
        friction = np.asarray(contact['tangential_force_by_feature_N'])
        force_by_feature = np.asarray(contact['force_by_feature_N'])
        active = np.asarray(contact['active'], bool)
        gaps = np.asarray(contact['gaps_m'])
        shoulder_index = contact['features'].index('SHOULDER')
        contact_radius = float(state['parameters']['dx']) / 2
        assert np.array_equal(active, gaps < contact_radius)
        assert np.array_equal(total, internal + contact_force)
        assert np.array_equal(contact_force, np.asarray(contact['force_N']))
        assert np.array_equal(force_by_feature, normal + friction)
        assert np.array_equal(contact_force, force_by_feature.sum(axis=1))
        assert np.array_equal(reaction, -total * (~free))
        residual = float(np.max(np.abs(total[free])))
        assert residual <= evidence['metrics']['force_limit_N']
        feature_rows = []
        for k, name in enumerate(contact['features']):
            feature_rows.append(dict(feature=name, active_particles=int(active[:, k].sum()),
                normal_resultant_on_material_N=normal[:, k].sum(axis=0).tolist(),
                friction_resultant_on_material_N=friction[:, k].sum(axis=0).tolist(),
                sum_normal_magnitudes_N=float(np.linalg.norm(normal[:, k], axis=1).sum())))
        crossed = np.isin(ids, state['bore_observation']['crossing_ids'])
        nearest_shoulder = int(np.argmin(gaps[:, shoulder_index]))
        crossed_indices = np.flatnonzero(crossed)
        nearest_crossed_shoulder = int(crossed_indices[np.argmin(gaps[crossed_indices, shoulder_index])])
        rr = np.linalg.norm(pos[:, :2], axis=1)
        radial = np.zeros_like(pos)
        radial[:, :2] = np.divide(pos[:, :2], rr[:, None], out=np.zeros_like(pos[:, :2]), where=rr[:, None] > 0)
        groups = {}
        bar = np.asarray(pa.component_id) == 1
        for name, mask in [('bar_above_tool_current', bar & (pos[:, 2] > .0025)),
                           ('bar_inside_bore_current', bar & (pos[:, 2] > .0015) & (pos[:, 2] < .0025)),
                           ('bar_under_bore_current', bar & (pos[:, 2] < .0015)),
                           ('crossed', crossed)]:
            groups[name] = dict(count=int(mask.sum()),
                contact_resultant_on_material_N=contact_force[mask].sum(axis=0).tolist(),
                sum_contact_radial_components_N=float(np.sum(contact_force[mask] * radial[mask])),
                sum_contact_magnitudes_N=float(np.linalg.norm(contact_force[mask], axis=1).sum()),
                internal_resultant_N=internal[mask].sum(axis=0).tolist())
        cuts = []
        for z_mm in [0, .5, 1, 1.5, 2, 2.5, 3, 4, 5]:
            upper = pos[:, 2] > z_mm / 1000
            cuts.append(dict(current_horizontal_plane_z_mm=z_mm, particles_above=int(upper.sum()),
                internal_resultant_on_upper_set_N=internal[upper].sum(axis=0).tolist(),
                contact_resultant_on_upper_set_N=contact_force[upper].sum(axis=0).tolist(),
                prescribed_reaction_on_upper_set_N=reaction[upper].sum(axis=0).tolist()))
        rows.append(dict(time_s=float(state['time_s']), step=int(state['step']),
            checkpoint=source, candidate_json=dict(path=str(jp.relative_to(ROOT)), sha256=sha(jp)),
            candidate_npz=dict(path=str(npz.relative_to(ROOT)), sha256=sha(npz)),
            interval_dt_s=motion['dt_s'], free_residual_N=residual,
            force_limit_N=evidence['metrics']['force_limit_N'],
            actuator_external_reaction_N=reaction[actuator].sum(axis=0).tolist(),
            nonactuator_external_reaction_N=reaction[~actuator].sum(axis=0).tolist(),
            tool_reaction_from_material_N=contact['particle_reaction_on_tool_N'],
            tool_reaction_torque_from_material_Nm=contact['particle_reaction_torque_on_tool_Nm'],
            prescribed_work_J=energy['prescribed_boundary_work_J'],
            tool_work_J=energy['tool_work_J'], features=feature_rows, groups=groups,
            horizontal_cuts=cuts,
            shoulder_geometry=dict(contact_radius_m=contact_radius,
                minimum_gap_all_particles_m=float(gaps[nearest_shoulder, shoulder_index]),
                minimum_clearance_to_contact_all_particles_m=float(gaps[nearest_shoulder, shoulder_index]-contact_radius),
                nearest_particle_id=int(ids[nearest_shoulder]),
                nearest_particle_component_id=int(pa.component_id[nearest_shoulder]),
                nearest_particle_initial_z_m=float(pa.z0[nearest_shoulder]),
                nearest_particle_current_position_m=pos[nearest_shoulder].tolist(),
                nearest_particle_is_crossed=bool(crossed[nearest_shoulder]),
                minimum_gap_crossed_particles_m=float(gaps[nearest_crossed_shoulder, shoulder_index]),
                minimum_clearance_to_contact_crossed_particles_m=float(gaps[nearest_crossed_shoulder, shoulder_index]-contact_radius),
                nearest_crossed_particle_id=int(ids[nearest_crossed_shoulder])),
            crossed_particles=[dict(particle_id=int(ids[k]),
                contact_force_N=contact_force[k].tolist(), internal_force_N=internal[k].tolist(),
                active_features=[contact['features'][j] for j in np.where(active[k])[0]],
                shoulder_gap_m=float(gaps[k, shoulder_index]),
                shoulder_clearance_to_contact_m=float(gaps[k, shoulder_index]-contact_radius),
                nearest_contact_feature=contact['features'][int(np.argmin(gaps[k]))],
                nearest_contact_feature_gap_m=float(np.min(gaps[k])),
                velocity_m_s=velocity[k].tolist()) for k in np.where(crossed)[0]]))
    output = dict(classification='READ_ONLY_SAVED_FINITE_STEP_CONTACT_EVIDENCE', physical_steps_executed=0,
        conventions=dict(positive_z='UP', force='Saved finite-step discrete normal and BE tangential forces on material, world coordinates.',
            reaction='External constraint force on material: -total_force*(~free). Tool reaction is opposite contact on material.',
            stress='Saved stress_Pa is deviatoric S; total Cauchy stress is S-K*(rho/rho0-1)*I.',
            residual='Quasi-static net free force near zero does not imply zero stress or absence of neighbor load.',
            cuts='Internal resultant summed over all material above a current horizontal plane. Pair antisymmetry cancels interactions wholly inside that set, leaving net interaction across the partition; this is a discrete resultant, not a calibrated continuum traction.',
            shoulder_geometry='gaps_m is distance from particle center to finite shoulder patch; contact is active when gap < radius=dx/2. clearance_to_contact=gap-radius. Nearest ID may change; these are geometric minima, not fixed-particle trajectories.',
            sampling='Parent plus 83 completed macro endpoints; selected final fine_2 interval at each endpoint, not all intermediate substeps.'),
        records=rows,
        summary=dict(sampled_intervals=len(rows),
            shoulder_max_active=max(x['features'][0]['active_particles'] for x in rows),
            shoulder_max_sum_normal_magnitudes_N=max(x['features'][0]['sum_normal_magnitudes_N'] for x in rows),
            all_particles_minimum_shoulder_gap_range_m=[min(x['shoulder_geometry']['minimum_gap_all_particles_m'] for x in rows),max(x['shoulder_geometry']['minimum_gap_all_particles_m'] for x in rows)],
            all_particles_nearest_shoulder_ids=sorted(set(x['shoulder_geometry']['nearest_particle_id'] for x in rows)),
            crossed_minimum_shoulder_gap_range_m=[min(x['shoulder_geometry']['minimum_gap_crossed_particles_m'] for x in rows),max(x['shoulder_geometry']['minimum_gap_crossed_particles_m'] for x in rows)],
            crossed_max_contact_force_norm_N=max(np.linalg.norm(p['contact_force_N']) for x in rows for p in x['crossed_particles']),
            actuator_z_reaction_range_N=[min(x['actuator_external_reaction_N'][2] for x in rows), max(x['actuator_external_reaction_N'][2] for x in rows)],
            bore_contact_z_resultant_range_N=[min(x['groups']['bar_inside_bore_current']['contact_resultant_on_material_N'][2] for x in rows), max(x['groups']['bar_inside_bore_current']['contact_resultant_on_material_N'][2] for x in rows)]),
        limitations=['Sampled intervals do not certify zero contact at every intervening substep.',
            'Force direction and kinematics do not prove a unique cause of limited outflow.',
            'No rerun, force reevaluation, resolution comparison, material validation or change to the frozen model.'])
    destination = Path(__file__).with_name('contact_readback.json')
    destination.write_text(json.dumps(output, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(output['summary'], indent=2))


if __name__ == '__main__':
    main()
