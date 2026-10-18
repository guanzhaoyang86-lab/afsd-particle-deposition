"""Named snapshot interface for the joint BE-thermal process-root experiment.

Mechanical accumulation reuses the old tested commit. New state classification
is explicit and cannot be read by the old continuation loader. File roundtrip
is not a two-trajectory recomputed restart validation.
"""
import hashlib
import numpy as np
from pysph.solver.utils import load
from .continuation_state import (commit as mechanical_commit, state as common_state,
                                 VALID_PROPERTIES)
from .updated_trial import particle_hash
from .process_thermal_trial import CLASSIFICATION as TRIAL_CLASSIFICATION
from .mechanics import positions, velocities, pressure, bulk_modulus
from .material import equivalent
from .output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CLASSIFICATION='STRUCTURAL_ROOT_BE_THERMAL_CONTINUATION_APPROXIMATE_NOT_PROCESS_VALIDATED'


def commit(pa,p,candidate,evaluation,checks):
    if evaluation.get('classification')!=TRIAL_CLASSIFICATION or checks.get('variant')!=TRIAL_CLASSIFICATION:
        raise ValueError('named joint process acceptance required')
    return mechanical_commit(pa,p,candidate,evaluation,checks)


def state(pa,p,time,stage_index,sources,config_hash,parent_reference,transport_cumulative,ledger):
    s=common_state(pa,p,time,stage_index,sources,config_hash,parent_reference,transport_cumulative,ledger)
    s.update(classification=CLASSIFICATION,trial_classification=TRIAL_CLASSIFICATION,
        temporal_accuracy='NOT_VALIDATED_LARGER_STEP_ROOT_DIAGNOSTIC_ONLY',
        ledger_scope='THIS_NAMED_VARIANT_ONLY; PARTICLE_PLASTIC_AND_RETURN_INCLUDE_PARENT',
        restart_scope='INTERNAL_ROUNDTRIP_ONLY; PUBLIC_RESUME_REQUIRES_PUBLISHED_HISTORY; RECOMPUTED_EQUIVALENCE_NOT_RUN')
    return s


def restore(path,p,sources,config_hash,expected_file_hash):
    """Internal roundtrip only; a lone NPZ does not authorize public resumption."""
    if hashlib.sha256(path.read_bytes()).hexdigest()!=expected_file_hash:raise ValueError('snapshot hash mismatch')
    data=load(str(path));pa=data['arrays']['solid'];s=data['solver_data']
    if (s.get('classification')!=CLASSIFICATION or s.get('trial_classification')!=TRIAL_CLASSIFICATION
            or s.get('accepted') is not True or s.get('restart_allowed') is not True
            or s.get('production_authorized') is not False):raise ValueError('not an accepted named snapshot')
    if s['parameters']!=p or s['source_sha256']!=sources or s['configuration_sha256']!=config_hash:
        raise ValueError('snapshot source/config mismatch')
    if s['particle_state_sha256']!=particle_hash(pa) or s['particle_ids']!=pa.particle_id.tolist():
        raise ValueError('snapshot ID/history mismatch')
    if len(np.unique(pa.particle_id))!=len(pa.x) or set(s['valid_properties'])!=set(VALID_PROPERTIES):
        raise ValueError('snapshot field contract mismatch')
    return pa,s


def save(directory,pa,p,s):
    if s.get('classification')!=CLASSIFICATION or s.get('particle_state_sha256')!=particle_hash(pa):
        raise ValueError('invalid named state for saving')
    checkpoint(directory/'accepted.npz',pa,s)
    S=pa.S.reshape(-1,3,3);pr=pressure(pa.rho,p['rho0'],bulk_modulus(p['shear'],p['poisson']))
    sigma=S-pr[:,None,None]*np.eye(3)
    scalars={k:getattr(pa,k).tolist() for k in ('particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','as_rs_tag')}
    scalars.update(mass_kg=pa.m.tolist(),density_kg_per_m3=pa.rho.tolist(),temperature_K=pa.temperature.tolist(),
        pressure_Pa=pr.tolist(),equivalent_stress_Pa=equivalent(S).tolist(),equivalent_plastic_strain=pa.ep.tolist(),
        cumulative_plastic_work_J=pa.plastic_work.tolist(),cumulative_return_diagnostic_J=pa.return_dissipation.tolist())
    for a,b,label in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        scalars['stress_'+label+'_Pa']=sigma[:,a,b].tolist()
    write_point_cloud(directory/'accepted.vtk',positions(pa).tolist(),scalars,
        dict(velocity_m_per_s=velocities(pa).tolist(),initial_position_m=np.column_stack((pa.x0,pa.y0,pa.z0)).tolist()),
        title=CLASSIFICATION+' time_s='+str(s['t']))
