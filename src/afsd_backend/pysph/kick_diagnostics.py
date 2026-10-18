"""Read-only observer of legacy candidate kicks; never supplies physics values."""
from collections import deque
from copy import deepcopy
from pathlib import Path
import json
import numpy as np
from .geometry import surfaces
from .mechanics import positions, velocities
from .output import checkpoint

FEATURES = ('BASE', 'SHOULDER', 'BORE_X_MINUS', 'BORE_X_PLUS', 'BORE_Y_MINUS', 'BORE_Y_PLUS')


def clean(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def active_keys(x, time, p, ids):
    return {(int(ids[idx]), feature) for feature, (indices, gaps, *_rest) in
            zip(FEATURES, surfaces(x, time, p)) for idx, gap in zip(indices, gaps)
            if gap < p['dx']/2}


class KickRecorder:
    def __init__(self, directory, provenance):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=False)
        self.provenance = provenance
        self.kicks = deque(maxlen=2)
        self.first_negative = None
        self.candidate = None
        self.accepted_pa = None
        self.accepted_state = None
        self.file = (self.directory/'contacts.jsonl').open('w', encoding='utf-8')
        write_json(self.directory/'provenance.json', provenance)

    def accept(self, solver):
        pa = solver.pa
        self.accepted_pa = pa.extract_particles(np.arange(pa.get_number_of_particles(), dtype=np.uint32))
        self.accepted_state = deepcopy(solver.state())
        self.accepted_state.update(classification='ACCEPTED_DEBUG', provenance=self.provenance)
        if self.kicks:
            self.kicks[-1]['status'] = 'ACCEPTED'
            self.file.write(json.dumps(clean(self.kicks[-1]), allow_nan=False)+'\n')
            self.file.flush()

    def record(self, solver, time, dt, internal_force, contact_force, old_v,
               acceleration, gradient, new_stress, new_ep, drho):
        pa, p = solver.pa, solver.p
        x = positions(pa)
        vf = old_v+dt*acceleration
        jint = dt*internal_force
        jcontact = dt*contact_force
        jconstraint = pa.m[:, None]*(vf-old_v)-jint-jcontact
        xend = x+.5*dt*(vf+np.column_stack((pa.ax, pa.ay, pa.az)))
        before = active_keys(positions(self.accepted_pa), solver.time, p, pa.particle_id)
        after = active_keys(xend, time+dt/2, p, pa.particle_id)
        records = []
        for feature, (indices, gaps, normals, wall_v, wall_t) in zip(FEATURES, surfaces(x, time, p)):
            for a, idx in enumerate(indices):
                f, _power, _potential, penetration = solver.contact.force(
                    gaps[a], normals[a], old_v[idx]-wall_v[a], p['dx']/2, pa.m[idx], dt)
                if penetration <= 0: continue
                n = normals[a]
                jn = dt*solver.kn*penetration*n
                jt = dt*f-jn
                point = x[idx]-gaps[a]*n
                center = np.array([p['traverse']*time, 0., p['layer']])
                omega = np.array([0., 0., 0. if feature == 'BASE' else p['omega']])
                translation = np.array([0. if feature == 'BASE' else p['traverse'], 0., 0.])
                wall_point_v = translation+np.cross(omega, point-center)
                material_point_v = old_v[idx]+gradient[idx]@(point-x[idx])
                vmid = .5*(old_v[idx]+vf[idx])
                tangent = np.eye(3)-np.outer(n, n)
                u0 = tangent@(old_v[idx]-wall_v[a])
                u1 = tangent@(vf[idx]-wall_v[a])
                umid = tangent@(vmid-wall_v[a])
                q = -float(jt@umid)
                pid = int(pa.particle_id[idx])
                m = pa.m[idx]
                v_int = old_v[idx]+jint[idx]/m
                v_contact = v_int+jcontact[idx]/m
                impulse_parts = dict(internal=jint[idx], normal=jn, friction=jt,
                    other_contacts=jcontact[idx]-jn-jt, constraint=jconstraint[idx])
                direction = u0/np.linalg.norm(u0) if np.linalg.norm(u0) > 0 else np.zeros(3)
                projections = {k: float(direction@(tangent@val)/m) for k, val in impulse_parts.items()}
                rounding = 64*np.finfo(float).eps*np.linalg.norm(jt)*(np.linalg.norm(old_v[idx])+np.linalg.norm(vf[idx])+2*np.linalg.norm(wall_v[a]))
                dk = .5*m*float(vf[idx]@vf[idx]-old_v[idx]@old_v[idx])
                work_internal = float(jint[idx]@vmid)
                work_contacts = float(jcontact[idx]@vmid)
                work_constraint = float(jconstraint[idx]@vmid)
                row = dict(step=solver.count+1, stage='ACCELERATION_CALLBACK_AFTER_CONSTRAINT_BEFORE_HEAT_CHECK',
                    particle_id=pid, contact_id=str(pid)+':'+feature, feature_id=feature, face_id=feature,
                    time_mid_s=time, dt_s=dt, position_mid_m=x[idx], contact_point=point,
                    force_application_point=x[idx], normal=n, gap_m=gaps[a], penetration_m=penetration,
                    active_before=(pid,feature) in before, active_mid=True, active_after=(pid,feature) in after,
                    geometry_after_definition='REJECTED_CANDIDATE_END_DRIFT_PREDICTION, not accepted geometry',
                    feature_frozen_during_kick=True,
                    velocities_definition='virtual partial impulse sums; actual integrator applies one simultaneous kick',
                    v_before_internal=old_v[idx], v_after_internal=v_int,
                    v_before_contact=v_int, v_after_contact=v_contact,
                    v_before_constraint=v_contact, v_after_constraint=vf[idx],
                    contact_force_evaluation_velocity=old_v[idx],
                    wall_translation_velocity=translation, omega=omega,
                    contact_point_wall_velocity=wall_point_v,
                    material_contact_point_velocity=material_point_v,
                    material_contact_point_velocity_status='AFFINE_RECONSTRUCTION_FOR_DIAGNOSTICS_ONLY',
                    actual_wall_velocity_used=wall_v[a], actual_material_velocity_used=old_v[idx],
                    relative_velocity=old_v[idx]-wall_v[a], tangential_relative_velocity=u0,
                    tangential_relative_midpoint=umid, tangential_relative_after=u1,
                    internal_impulse=jint[idx], normal_contact_impulse=jn, tangential_contact_impulse=jt,
                    constraint_impulse=jconstraint[idx], actuator=bool(pa.actuator[idx]),
                    impulse_velocity_projection_on_initial_slip=projections,
                    initial_slip_speed=float(np.linalg.norm(u0)),
                    endpoint_projection=float(direction@u1), midpoint_projection=float(direction@umid),
                    friction_work_raw=q, friction_work_definition='positive dissipation -Jt dot relative_midpoint, joules',
                    friction_lab_work=float(jt@vmid),
                    normal_contact_work=float(jn@vmid), constraint_work=work_constraint,
                    constraint_sequential_virtual_work=float(jconstraint[idx]@(.5*(v_contact+vf[idx]))),
                    internal_work=work_internal, all_contacts_lab_work=work_contacts,
                    delta_kinetic_energy=dk,
                    local_energy_residual=dk-work_internal-work_contacts-work_constraint,
                    local_energy_residual_definition='kinetic kick work identity only, NOT material free energy',
                    floating_point_work_scale_J=rounding,
                    negative_vs_roundoff=abs(q)/rounding if rounding > 0 else None,
                    heat_guard_trigger=q < -1e-12*dt)
                records.append(row)
        for row in records:
            others = [r for r in records if r['particle_id']==row['particle_id'] and r['contact_id']!=row['contact_id']]
            row['impulse_from_each_other_contact'] = [dict(contact_id=r['contact_id'],
                normal_impulse=r['normal_contact_impulse'], tangent_impulse=r['tangential_contact_impulse']) for r in others]
            row['features_before'] = sorted(f for pid,f in before if pid == row['particle_id'])
            row['features_after'] = sorted(f for pid,f in after if pid == row['particle_id'])
            if row['heat_guard_trigger'] and self.first_negative is None:
                self.first_negative = deepcopy(row)
        kick = dict(step=solver.count+1, status='CANDIDATE', records=records, active_contacts=len(records))
        self.kicks.append(kick)
        self.candidate = dict(position_mid=x.copy(), position_end_predicted=xend,
            velocity_before=old_v.copy(), velocity_after_predicted=vf,
            stress_candidate=new_stress.copy(), ep_candidate=new_ep.copy(),
            rho_candidate=pa.rho+dt*drho, temperature_not_advanced=pa.temperature.copy(),
            particle_id=pa.particle_id.copy(), time_mid_s=time, time_end_s=time+dt/2)

    def reject(self, solver, error):
        if self.kicks:
            self.kicks[-1]['status'] = 'REJECTED'
            self.file.write(json.dumps(clean(self.kicks[-1]), allow_nan=False)+'\n')
        self.file.close()
        write_json(self.directory/'last_two_kicks.json', list(self.kicks))
        write_json(self.directory/'first_negative_event.json', self.first_negative)
        checkpoint(self.directory/'last_accepted.npz', self.accepted_pa, self.accepted_state)
        checkpoint(self.directory/'REJECTED_partial_integrator_state.npz', solver.pa,
            {**deepcopy(solver.state()), 'classification':'REJECTED', 'restart_allowed':False,
             'reason':str(error), 'provenance':self.provenance, 'stage':'after_drift_before_actual_kick'})
        if self.candidate:
            np.savez_compressed(self.directory/'REJECTED_candidate_fields.npz', **self.candidate)
        write_json(self.directory/'REJECTED_manifest.json', dict(classification='REJECTED',
            restart_allowed=False, not_valid_physical_endpoint=True, error=str(error),
            candidate_fields_status='predicted kick fields, thermal NOT_EVALUATED; no synthetic final field',
            last_accepted_step=self.accepted_state['count'], first_negative=self.first_negative,
            provenance=self.provenance))
