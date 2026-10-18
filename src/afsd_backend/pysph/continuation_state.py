"""Accepted updated-stage history and strict same-source snapshot loading.

No force/material/thermal equations are implemented here. Scratch fields are
not presented as computed internal energy; real increments remain separate.
"""
import copy
import hashlib
import numpy as np
from pysph.solver.utils import load
from .updated_trial import particle_hash
from .mechanics import positions, velocities, pressure, bulk_modulus
from .material import equivalent
from .geometry import tool_pose
from .output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CLASSIFICATION='UPDATED_STAGE_ANDERSON_CONTINUATION_APPROXIMATE'
VALID_PROPERTIES=('x','y','z','u','v','w','rho','temperature','S','D','ep','ep_rate','strain_rate',
    'plastic_work','return_dissipation','particle_id','material_id','component_id','tracer_type',
    'wholebar_tag','actuator','as_rs_tag','m','h','x0','y0','z0','feed_crossed','p')


def commit(pa, p, candidate, evaluation, checks):
    """Copy a validated candidate, accumulating real accepted increments once."""
    if checks.get('status')!='ACCEPTED_NUMERICAL_STAGE_APPROXIMATE_NOT_ENERGY_CERTIFIED':
        raise ValueError('complete numerical stage acceptance required')
    if evaluation['accepted_state_hash']!=particle_hash(pa): raise ValueError('commit initial history mismatch')
    if not np.array_equal(pa.particle_id,evaluation['particle_id']): raise ValueError('commit ID mismatch')
    out=pa.extract_particles(np.arange(len(pa.x),dtype=np.uint32)); dt=evaluation['dt_s']
    out.x[:],out.y[:],out.z[:]=np.asarray(candidate['position_m']).T
    out.u[:],out.v[:],out.w[:]=np.asarray(candidate['velocity_m_s']).T
    out.rho[:]=candidate['density_kg_m3']; out.temperature[:]=candidate['temperature_K']
    material=evaluation['material']; out.S[:]=np.asarray(material['deviatoric_stress_Pa']).ravel()
    out.ep[:]=material['ep']; out.ep_rate[:]=(out.ep-pa.ep)/dt; out.strain_rate[:]=material['strain_rate_per_s']
    gradient=np.asarray(evaluation['midpoint']['velocity_gradient_per_s'])
    out.D[:]=(.5*(gradient+gradient.swapaxes(1,2))).ravel()
    out.plastic_work[:]+=material['plastic_work_J']; out.return_dissipation[:]+=material['return_dissipation_J']
    out.au[:],out.av[:],out.aw[:]=((velocities(out)-velocities(pa))/dt).T
    out.arho[:]=(out.rho-pa.rho)/dt
    out.p[:]=pressure(out.rho,p['rho0'],bulk_modulus(p['shear'],p['poisson']))
    crossed=(out.component_id==1)&(pa.z>=p['layer'])&(out.z<p['layer'])
    out.feed_crossed[crossed]=1.
    for key in VALID_PROPERTIES:
        if not np.isfinite(getattr(out,key)).all(): raise ValueError('nonfinite accepted property '+key)
    return out


def state(pa, p, time, stage_index, sources, config_hash, parent_reference, transport_cumulative, ledger):
    centre, angle=tool_pose(time,p['traverse'],p['omega'],p['layer'])
    return dict(classification=CLASSIFICATION,accepted=True,restart_allowed=True,production_authorized=False,
        t=float(time),stage_index=int(stage_index),parameters=copy.deepcopy(p),source_sha256=dict(sources),
        configuration_sha256=config_hash,particle_state_sha256=particle_hash(pa),
        particle_ids=pa.particle_id.tolist(),parent_reference=copy.deepcopy(parent_reference),
        transport_cumulative=copy.deepcopy(transport_cumulative),ledger=copy.deepcopy(ledger),
        ledger_scope='CONTINUATION_ONLY; PARTICLE_CUMULATIVE_PLASTIC_AND_RETURN_INCLUDE_PARENT_HISTORY',
        tool_pose=dict(centre_m=centre.tolist(),angle_rad=float(angle)),valid_properties=list(VALID_PROPERTIES),
        scratch_not_solved=['e','ae','ax','ay','az'],
        restart_scope='STRICT_NEW_LOADER_ONLY; REAL_RECOMPUTED_RESTART_EQUIVALENCE_NOT_YET_CERTIFIED',
        total_free_energy_status='NOT_CERTIFIED_DIAGNOSTIC_ONLY')


def restore(path,p,sources,config_hash,expected_file_hash):
    """Internal exact roundtrip; public resumption additionally needs a published history entry."""
    if hashlib.sha256(path.read_bytes()).hexdigest()!=expected_file_hash: raise ValueError('snapshot file hash mismatch')
    data=load(str(path)); pa=data['arrays']['solid']; s=data['solver_data']
    if (s.get('classification')!=CLASSIFICATION or s.get('accepted') is not True or
            s.get('restart_allowed') is not True or s.get('production_authorized') is not False):
        raise ValueError('not an accepted continuation snapshot')
    if s['parameters']!=p or s['source_sha256']!=sources or s['configuration_sha256']!=config_hash:
        raise ValueError('snapshot parameter/source/configuration mismatch')
    if s['particle_state_sha256']!=particle_hash(pa) or s['particle_ids']!=pa.particle_id.tolist():
        raise ValueError('snapshot persistent history/ID mismatch')
    if len(np.unique(pa.particle_id))!=len(pa.x): raise ValueError('duplicate snapshot IDs')
    if set(s['valid_properties'])!=set(VALID_PROPERTIES): raise ValueError('snapshot field contract mismatch')
    return pa,s


def save(directory,pa,p,s):
    checkpoint(directory/'accepted.npz',pa,s)
    S=pa.S.reshape(-1,3,3); pr=pressure(pa.rho,p['rho0'],bulk_modulus(p['shear'],p['poisson']))
    sigma=S-pr[:,None,None]*np.eye(3)
    scalars={k:getattr(pa,k).tolist() for k in ('particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','as_rs_tag')}
    scalars.update(mass_kg=pa.m.tolist(),density_kg_per_m3=pa.rho.tolist(),temperature_K=pa.temperature.tolist(),
        pressure_Pa=pr.tolist(),equivalent_stress_Pa=equivalent(S).tolist(),equivalent_plastic_strain=pa.ep.tolist(),
        cumulative_plastic_work_J=pa.plastic_work.tolist(),cumulative_return_diagnostic_J=pa.return_dissipation.tolist())
    for a,b,label in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        scalars['stress_'+label+'_Pa']=sigma[:,a,b].tolist()
    write_point_cloud(directory/'accepted.vtk',positions(pa).tolist(),scalars,
        dict(velocity_m_per_s=velocities(pa).tolist(),initial_position_m=np.column_stack((pa.x0,pa.y0,pa.z0)).tolist()),
        title=CLASSIFICATION+' time_s='+str(s['t'])+'; NOT_DEPOSITION_VALIDATION')
