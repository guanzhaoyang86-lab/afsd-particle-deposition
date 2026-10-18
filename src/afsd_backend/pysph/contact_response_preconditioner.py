"""Assemble a named fixed-node NUMERICAL response from an audited candidate.

This is not a continuum cross-zero tangent or a new contact force. All actual
candidate acceptance remains with the unmodified nonlinear residual/guards.
"""
import numpy as np
from .contact_path_response import path_response

STRATEGY='FIXED_QUADRATURE_DISCRETE_CONTACT_PRECONDITIONER_EXPERIMENT_V1'


def assemble(contact,right_boundary,endpoint_velocity,dt,p):
    velocity=np.asarray(endpoint_velocity,float)
    if velocity.ndim!=2 or velocity.shape[1]!=3 or not np.isfinite(velocity).all():
        raise ValueError('finite Cartesian endpoint velocities required')
    count=len(velocity);matrix=np.zeros((3*count,3*count));details=[]
    right={row['contact_id']:row for row in right_boundary['contacts']}
    if len(right)!=len(right_boundary['contacts']):raise ValueError('duplicate contact identity')
    for row in contact['contacts']:
        cid=row['contact_id'];end=right[cid];index=row['particle_index'];normal=np.asarray(row['normal'])
        if not 0<=index<count or end['particle_index']!=index or end['particle_id']!=row['particle_id']:
            raise ValueError('contact/particle identity mismatch')
        if not np.array_equal(normal,end['normal']):raise ValueError('changed normal unsupported')
        certificate=row.get('event_path_certificate')
        mode=None if certificate is None else certificate['mode']
        if mode not in (None,'private_continuation'):
            raise ValueError('zero endpoint derivative outside registered finite-node experiment')
        kn=end['normal_stiffness_N_m'];mu=p['friction']
        def make(fn,u,force):
            speed=float(np.linalg.norm(u));sign=1
            if mode=='private_continuation':
                if speed<=certificate['zero_slip_bound_m_s']:
                    raise ValueError('zero-override node has no registered smooth response')
                sign=1 if np.dot(u,certificate['incoming_direction'])>=0 else -1
            if speed==0:raise ValueError('zero slip has no pointwise response')
            expected=-sign*mu*fn*u/speed
            bound=512*np.finfo(float).eps*max(mu*fn,np.linalg.norm(force),np.finfo(float).tiny)
            if np.linalg.norm(expected-force)>bound:raise ValueError('response branch differs from actual force')
            return dict(normal=normal,tangential_relative_velocity_m_s=u,normal_force_N=fn,
                normal_stiffness_N_m=kn,friction=mu,branch_sign=sign)
        nodes=[]
        for node in row['quadrature']:
            item=make(node['normal_force_N'],np.asarray(node['tangential_relative_velocity_m_s']),
                np.asarray(node['tangential_force_N']))
            item.update(fraction=node['fraction'],weight=node['weight']);nodes.append(item)
        relative=velocity[index]-end['wall_velocity_m_s'];slip=relative-(relative@normal)*normal
        endpoint=make(end['normal_force_N'],slip,end['tangential_force_N'])
        response=path_response(nodes,endpoint,dt,p.get('hht_weight',.9))
        block=response['root_contact_response_kg'];sl=slice(3*index,3*index+3)
        matrix[sl,sl]+=block
        eigen=np.linalg.eigvals(block)
        details.append(dict(contact_id=cid,particle_id=row['particle_id'],particle_index=index,
            quadrature_order=len(nodes),response=response,block_eigenvalues_real_kg=eigen.real,
            block_eigenvalues_imaginary_kg=eigen.imag,mode=mode))
    return matrix,dict(strategy=STRATEGY,contacts=details,
        matrix_norm_kg=float(np.linalg.norm(matrix)),matrix_asymmetry_norm_kg=float(np.linalg.norm(matrix-matrix.T)),
        accepted=False,restart_allowed=False,can_publish=False,only_search_direction_changed=True,
        continuous_cross_zero_jacobian_certified=False,adaptive_quadrature_derivative_included=False,
        old_physical_and_root_gates_unchanged=True)
