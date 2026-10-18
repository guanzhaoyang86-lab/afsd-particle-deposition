"""Independent array/journal accounts; NEVER candidate acceptance or physics."""
import math
import numpy as np


def _json(value):
    if isinstance(value, dict): return {str(k):_json(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [_json(v) for v in value]
    if isinstance(value, np.ndarray): return _json(value.tolist())
    if isinstance(value, np.generic): return _json(value.item())
    if isinstance(value, float) and not math.isfinite(value): return str(value)
    return value


class _Checks:
    def __init__(self): self.rows={}

    def compare(self,name,actual,reference,atol=0.,rtol=0.,context=None):
        a=np.asarray(actual,float);b=np.asarray(reference,float)
        if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError('invalid/mismatched/nonfinite '+name)
        limits=np.broadcast_to(np.asarray(atol,float),a.shape)+rtol*np.abs(b)
        if not np.isfinite(limits).all() or np.any(limits<0):raise ValueError('invalid bound '+name)
        error=np.abs(a-b);passed=bool(np.all(error<=limits))
        ratios=np.divide(error,limits,out=np.zeros_like(error),where=limits>0)
        ratios[(limits==0)&(error>0)]=np.inf
        self._record(name,passed,float(np.max(error,initial=0)),float(np.max(limits,initial=0)),
                     float(np.max(ratios,initial=0)),context)

    def require(self,name,condition,context=None):
        self._record(name,bool(condition),0. if condition else 1.,0.,0. if condition else np.inf,context)

    def _record(self,name,passed,error,bound,ratio,context):
        row=self.rows.setdefault(name,dict(name=name,count=0,passed=True,maximum_error=0.,
            maximum_bound=0.,maximum_error_over_bound=0.,first_failure=None))
        row['count']+=1;row['passed']=row['passed'] and passed
        row['maximum_error']=max(row['maximum_error'],error);row['maximum_bound']=max(row['maximum_bound'],bound)
        row['maximum_error_over_bound']=max(row['maximum_error_over_bound'],ratio)
        if not passed and row['first_failure'] is None:row['first_failure']=context or {}


def audit(pa,p,t0,dt,evaluation):
    """Return finite JSON data; inspect map accounts, never authorize a state.

The parent is a read-only PA-like object exposing m/rho/u/v/w/temperature/ep,
particle_id and x/y/z arrays. All field names otherwise follow the frozen
ProcessThermalTrial journal. No root convergence is required or implied.
"""
    checks=_Checks();summary={};out=dict(passed=False,status='FAIL_MAP_INTERNAL_ACCOUNTS',
        classification='INDEPENDENT_ARRAY_MAP_ACCOUNTS_NOT_CANDIDATE_ACCEPTANCE',
        accepted=False,restart_allowed=False,can_publish=False,production_authorized=False,
        checks=[],summary=summary)
    def array(value,shape,name):
        result=np.asarray(value,float)
        if result.shape!=shape or not np.isfinite(result).all():raise ValueError('invalid '+name)
        return result
    def work(name,actual,reference,scale,context=None):
        checks.compare(name,actual,reference,1e-15+1e-9*float(scale),context=context)
    try:
        if not np.isfinite([t0,dt]).all() or t0<0 or dt<=0:raise ValueError('positive finite time interval required')
        e=evaluation;m=np.asarray(pa.m,float);n=len(m);ids=np.asarray(pa.particle_id)
        if m.shape!=(n,) or n==0 or not np.isfinite(m).all() or np.any(m<=0):raise ValueError('positive mass vector required')
        v0=np.column_stack((pa.u,pa.v,pa.w));x0=np.column_stack((pa.x,pa.y,pa.z))
        v0=array(v0,(n,3),'parent velocity');temperature0=array(pa.temperature,(n,),'parent T')
        ep0=array(pa.ep,(n,),'parent ep');oldrho=array(pa.rho,(n,),'parent rho')
        checks.require('uncommitted_BE_metadata',e['classification']=='PROCESS_GLOBAL_EVENT_AVERAGE_BE_THERMAL_APPROXIMATE'
            and e['accepted'] is False and e['restart_allowed'] is False and e['production_authorized'] is False
            and e['temperature_layers']=='FKS_TMID; PAIR_AND_WALL_CONDUCTION_T1_BE; INCREMENT_FROM_T0')
        checks.require('particle_identity',np.array_equal(e['particle_id'],ids) and len(np.unique(ids))==n)
        checks.compare('time_before',e['time_before_s'],t0);checks.compare('dt',e['dt_s'],dt)
        checks.compare('time_after',e['time_after_s'],t0+dt);checks.compare('time_mid',e['time_mid_s'],t0+.5*dt)
        mapped=e['map'];mid=e['midpoint'];mat=e['material'];boundary=e['boundary'];thermal=e['thermal']
        vmap=array(mapped['velocity_m_s'],(n,3),'map velocity');vm=array(mapped['mean_velocity_m_s'],(n,3),'map mean')
        tmap=array(mapped['temperature_K'],(n,),'map T');force=array(e['internal_force_N'],(n,3),'internal force')
        checks.compare('map_boundary_velocity',vmap,boundary['velocity'])
        checks.compare('map_boundary_mean',vm,boundary['work_velocity'])
        checks.compare('map_position_identity',mapped['position_m'],x0+dt*vm,1e-15,1e-13)
        ep=array(mat['ep'],(n,),'material ep');stress=array(mat['deviatoric_stress_Pa'],(n,3,3),'material stress')
        rhomid=array(mid['density_kg_m3'],(n,),'midpoint rho')
        if np.any(rhomid<=0) or np.any(oldrho<=0):raise ValueError('positive physical density required')
        dep=ep-ep0;q=np.sqrt(1.5*np.sum(stress*stress,axis=(1,2)))
        plastic=m/rhomid*q*dep;returned=m/rhomid*1.5*p['shear']*dep*dep
        checks.require('plastic_history_monotone',np.all(dep>=0))
        checks.compare('plastic_work_identity',mat['plastic_work_J'],plastic,1e-18,1e-10)
        checks.compare('return_dissipation_separate',mat['return_dissipation_J'],returned,1e-18,1e-10)
        rate=array(mat['strain_rate_per_s'],(n,),'strain rate')
        checks.require('finite_nonnegative_rate',np.all(rate>=0))
        positive=np.flatnonzero(dep>0);summary.update(particle_count=n,plastic_particle_ids=ids[positive].tolist(),
            positive_plastic_particle_count=len(positive),plastic_work_J=float(math.fsum(plastic)),
            return_dissipation_J=float(math.fsum(returned)),positive_plastic_work_J=float(math.fsum(plastic[positive])),
            plastic_particles=[dict(particle_id=int(ids[k]),delta_ep=float(dep[k]),q_new_Pa=float(q[k]),
                plastic_work_J=float(plastic[k]),strain_rate_per_s=float(rate[k]),
                fks_rate_floor_active=bool(rate[k]<p['rate_floor'])) for k in positive],
            fks_rate_floor_particle_count=int(np.count_nonzero(rate<p['rate_floor'])))
        eps=np.finfo(float).eps;seen=set();wallrow=np.zeros(n);bath=np.zeros(n);qpart=np.zeros(n);qbound=np.zeros(n)
        friction_terms=np.zeros(n,dtype=np.int64);face_counts=np.zeros(n,dtype=np.int64)
        contact_ids=[];features=[];kinds={};negative=[];events=[];subintervals=0;multi=0;active_actuators=0;nonzero_reactions=0
        totals={name:0. for name in ['friction_raw_J','normal_material_work_J','normal_wall_work_J',
            'normal_relative_work_J','wall_work_J','actuator_work_J','internal_work_J','delta_kinetic_J']}
        maxmomentum=0.;maximum_free_reaction=0.;maximum_constraint_error=0.
        for trace in boundary['traces']:
            idx=int(trace['particle_index']);pid=int(trace['particle_id'])
            if not 0<=idx<n or idx in seen or pid!=ids[idx]:raise ValueError('duplicate/mismatched trace identity')
            seen.add(idx);ctx=dict(particle_id=pid,particle_index=idx);faces=trace['contact_rows'];nc=len(faces)
            face_counts[idx]=nc
            if nc>1:multi+=1
            active_actuators+=bool(trace['actuator_active'])
            checks.compare('trace_mass',trace['mass_kg'],m[idx],context=ctx)
            checks.compare('trace_initial_velocity',trace['material_velocity_before_m_s'],v0[idx],context=ctx)
            checks.compare('trace_internal_force',trace['internal_force_N'],force[idx],context=ctx)
            checks.compare('trace_dt',trace['dt_s'],dt,context=ctx)
            checks.compare('trace_mid_time',trace['time_mid_s'],t0+.5*dt,context=ctx)
            free=array(trace['free_axes'],(3,),'free axes');prescribed=array(trace['prescribed_velocity'],(3,),'prescribed velocity')
            checks.require('binary_free_axes',np.all((free==0)|(free==1)),ctx)
            normals=np.asarray([r['normal'] for r in faces],float).reshape(nc,3)
            walls=np.asarray([r['wall_velocity'] for r in faces],float).reshape(nc,3)
            fn=np.asarray([r['normal_force_N'] for r in faces],float)
            checks.compare('unit_normals',np.linalg.norm(normals,axis=1),np.ones(nc),1e-12,context=ctx)
            checks.require('active_contact_force',np.isfinite(fn).all() and np.all(fn>0),ctx)
            rowids=[str(f['contact_id']) for f in faces]
            checks.require('unique_contact_per_trace',len(set(rowids))==nc,ctx)
            contact_ids.extend(rowids)
            for face in faces:
                feature=int(face['feature_id']);features.append(feature)
                twall=p['base_temperature'] if feature==0 else p['tool_temperature']
                checks.compare('reservoir_temperature',face['wall_temperature_K'],twall,context=dict(ctx,contact_id=str(face['contact_id'])))
                wallrow[idx]+=p['contact_conductance']*p['dx']**2
                bath[idx]+=p['contact_conductance']*p['dx']**2*twall
            normalforce=fn[:,None]*normals;frame=walls[0] if nc else v0[idx]
            speed=max(np.linalg.norm(v0[idx]-frame),float(np.linalg.norm(walls-frame)) if nc else 0.,
                dt*np.linalg.norm(free*(force[idx]+normalforce.sum(axis=0))/m[idx]),np.finfo(float).tiny)
            result=trace['result'];covered=0.;previous=v0[idx];mean=np.zeros(3)
            jn_total=np.zeros((nc,3));jt_total=np.zeros((nc,3));jc_total=np.zeros(3);face_q=np.zeros(nc)
            traceworks={key:0. for key in ['other_work_J','constraint_work_J','wall_work_J','delta_kinetic_J']}
            for ordinal,interval in enumerate(result['intervals']):
                subintervals+=1;ic=dict(ctx,interval=ordinal,contact_ids=rowids)
                friction_terms[idx]+=nc
                h=float(interval['dt_s']);start=float(interval['start_s'])
                if not np.isfinite([h,start]).all() or h<=0:raise ValueError('invalid subinterval')
                checks.compare('subinterval_start',start,covered,128*eps*dt,context=ic)
                before=array(interval['v_before'],(3,),'interval before');after=array(interval['v_after'],(3,),'interval after')
                checks.compare('subinterval_history',before,previous,1e-15,1e-13,ic)
                midpoint=.5*(before+after);rel0=before-frame;rel1=after-frame;wr=walls-frame
                relative=.5*(rel0+rel1)-wr;tangent=relative-normals*np.sum(relative*normals,axis=1)[:,None]
                jn=array(np.asarray(interval['normal_impulse']).reshape(nc,3),(nc,3),'normal impulse')
                jt=array(np.asarray(interval['tangential_impulse']).reshape(nc,3),(nc,3),'tangent impulse')
                jc=array(interval['constraint_impulse'],(3,),'constraint impulse')
                raw=array(interval['friction_heat_raw_J'],(nc,),'raw face heat');bound=array(interval['friction_rounding_bound_J'],(nc,),'rounding bound')
                expected_bound=256*eps*np.linalg.norm(jt,axis=1)*(speed+np.linalg.norm(rel0)+np.linalg.norm(rel1)+2*np.linalg.norm(wr,axis=1))
                checks.compare('original_reference_frame',interval['reference_frame_velocity'],frame,context=ic)
                checks.compare('original_speed_scale',interval['relative_speed_scale_m_s'],speed,0.,1e-13,ic)
                checks.compare('original_rounding_certificate',bound,expected_bound,1e-30,1e-12,ic)
                checks.compare('tangent_midpoint',np.asarray(interval['tangential_relative_midpoint']).reshape(nc,3),tangent,
                    128*eps*(speed+np.linalg.norm(rel0)+np.linalg.norm(rel1)),context=ic)
                raw_ref=-np.sum(jt*tangent,axis=1)
                checks.compare('raw_friction_identity',raw,raw_ref,1e-25,1e-12,ic)
                checks.require('individual_raw_friction_not_clipped',np.all(bound>=0) and np.all(raw>=-bound),ic)
                for a in np.flatnonzero(raw<0):negative.append(dict(ic,contact_id=rowids[a],feature_id=int(faces[a]['feature_id']),
                    raw_Q_J=float(raw[a]),original_roundoff_J=float(bound[a])))
                momentum=m[idx]*(after-before)-h*force[idx]-jn.sum(axis=0)-jt.sum(axis=0)-jc
                checks.compare('interval_momentum',momentum,np.zeros(3),m[idx]*1e-9+1e-17,context=ic)
                checks.compare('normal_impulse_identity',jn,h*normalforce,m[idx]*1e-9+1e-17,context=ic)
                checks.compare('interval_free_reaction',free*jc,np.zeros(3),1e-17,context=ic)
                wc=float(jc@midpoint);wi=float(h*force[idx]@midpoint);wn=jn@midpoint
                ww=float(np.sum((jn+jt)*walls));dk=float(.5*m[idx]*((after-before)@(after+before)))
                scale=abs(wc)+abs(wi)+float(np.abs(wn).sum())+float(np.abs(jt@midpoint).sum())+abs(dk)
                for name,actual,reference in [('constraint_work',interval['constraint_work_J'],wc),
                    ('internal_work',interval['other_work_J'],wi),('normal_work',interval['normal_work_J'],wn),
                    ('wall_work',interval['wall_work_J'],ww),('delta_kinetic',interval['delta_kinetic_J'],dk)]:work(name,actual,reference,scale,ic)
                work('interval_kinetic_identity',dk,wc+wi+float(np.sum((jn+jt)*midpoint)),scale,ic)
                normalwall=float(np.sum(jn*walls));normalmaterial=float(wn.sum())
                work('interval_contact_work_identity',dk+float(raw.sum()),
                    wc+wi+ww+normalmaterial-normalwall,scale+abs(ww)+float(np.abs(raw).sum()),ic)
                for key,val in [('friction_raw_J',raw.sum()),('normal_material_work_J',normalmaterial),
                    ('normal_wall_work_J',normalwall),('normal_relative_work_J',normalmaterial-normalwall),
                    ('wall_work_J',ww),('actuator_work_J',wc),('internal_work_J',wi),('delta_kinetic_J',dk)]:totals[key]+=float(val)
                for key,val in [('other_work_J',wi),('constraint_work_J',wc),('wall_work_J',ww),('delta_kinetic_J',dk)]:traceworks[key]+=val
                kind=str(interval['kind']);kinds[kind]=kinds.get(kind,0)+1
                jn_total+=jn;jt_total+=jt;jc_total+=jc;face_q+=raw;qbound[idx]+=float(bound.sum())
                mean+=h*midpoint/dt;covered+=h;previous=after
            checks.compare('trace_coverage',covered,dt,128*eps*dt,context=ctx)
            checks.compare('trace_final_velocity',result['velocity'],previous,1e-15,1e-13,ctx)
            checks.compare('trace_map_velocity',vmap[idx],previous,1e-15,1e-13,ctx)
            checks.compare('trace_map_mean',vm[idx],mean,1e-9,context=ctx)
            checks.compare('trace_normal_impulse',np.asarray(result['normal_impulse']).reshape(nc,3),jn_total,m[idx]*1e-9+1e-17,context=ctx)
            checks.compare('trace_tangent_impulse',np.asarray(result['tangential_impulse']).reshape(nc,3),jt_total,m[idx]*1e-9+1e-17,context=ctx)
            checks.compare('trace_constraint_impulse',result['constraint_impulse'],jc_total,m[idx]*1e-9+1e-17,context=ctx)
            checks.compare('trace_friction_sum',result['friction_heat_J'],face_q,1e-25,1e-12,ctx)
            for key,value in traceworks.items():work('trace_'+key,result[key],value,sum(abs(v) for v in traceworks.values()),ctx)
            defect=m[idx]*(vmap[idx]-v0[idx])-dt*force[idx]-jn_total.sum(axis=0)-jt_total.sum(axis=0)-jc_total
            checks.compare('total_particle_momentum',defect,np.zeros(3),m[idx]*1e-9+1e-17,context=ctx)
            checks.compare('total_free_reaction',free*jc_total,np.zeros(3),1e-17,context=ctx)
            constrained=(1-free)*(vmap[idx]-prescribed)
            checks.compare('mapped_prescribed_velocity',constrained,np.zeros(3),1e-12,context=ctx)
            nonzero_reactions+=bool(np.linalg.norm((1-free)*jc_total)>0)
            maxmomentum=max(maxmomentum,float(np.max(abs(defect))));maximum_free_reaction=max(maximum_free_reaction,float(np.max(abs(free*jc_total))))
            maximum_constraint_error=max(maximum_constraint_error,float(np.max(abs(constrained))))
            qpart[idx]=float(face_q.sum())
            events.extend(dict(particle_id=pid,**_json(event)) for event in result['events'])
        checks.require('all_particle_traces',len(seen)==n)
        checks.compare('boundary_friction_once',boundary['friction_raw_J'],qpart,1e-25,1e-12)
        cert=e['thermal_boundary_certificate']
        for key,ref,atol,rtol in [('wall_row_W_K',wallrow,1e-15,1e-13),
            ('wall_weighted_temperature_W',bath,1e-15,1e-13),('friction_rounding_bound_J',qbound,1e-30,1e-12),
            ('friction_raw_J',qpart,1e-25,1e-12)]:checks.compare('thermal_certificate_'+key,cert[key],ref,atol,rtol)
        checks.compare('boundary_wall_row',boundary['wall_row_W_K'],wallrow,1e-15,1e-13)
        workscale=sum(abs(v) for v in totals.values())
        for name,value in totals.items():work('boundary_ledger_'+name,boundary['ledger'][name],value,workscale)
        # Independent radial undirected conductance, no thermal/kernel callback.
        first=np.asarray(mid['i']);second=np.asarray(mid['j']);npairs=len(first)
        if first.shape!=(npairs,) or second.shape!=first.shape or first.dtype.kind not in 'iu' or second.dtype.kind not in 'iu':raise ValueError('integer pair indices required')
        if np.any(first<0) or np.any(second<0) or np.any(first>=n) or np.any(second>=n) or np.any(first==second):raise ValueError('invalid pair indices')
        checks.require('unique_undirected_pairs',len(np.unique(np.sort(np.column_stack((first,second)),axis=1),axis=0))==npairs)
        xmid=array(mid['position_m'],(n,3),'midpoint position');g=array(mid['g'],(npairs,3),'raw pair gradient')
        rij=xmid[first]-xmid[second];r2=np.sum(rij*rij,axis=1)
        if np.any(r2<=0):raise ValueError('coincident thermal pair')
        conductivity=np.broadcast_to(np.asarray(p['conductivity'],float),(n,))
        if np.any(conductivity<=0) or not np.isfinite(conductivity).all():raise ValueError('positive conductivity required')
        conductance=m[first]*m[second]/(rhomid[first]*rhomid[second])*(4*conductivity[first]*conductivity[second]/(conductivity[first]+conductivity[second]))*(-np.sum(rij*g,axis=1))/r2
        checks.require('nonnegative_conductance',np.all(conductance>=0))
        checks.compare('pair_conductance_reconstruction',thermal['pair_conductance_W_K'],conductance,1e-25,1e-12)
        checks.compare('thermal_map_temperature',thermal['temperature'],tmap)
        pair=dt*conductance*(tmap[second]-tmap[first]);node=np.zeros(n)
        np.add.at(node,first,pair);np.add.at(node,second,-pair)
        wall=dt*(bath-wallrow*tmap);chi=float(p['plastic_heat_fraction']);partition=float(p['heat_partition'])
        cp=np.broadcast_to(np.asarray(p['cp'],float),(n,))
        if not 0<=chi<=1 or not 0<=partition<=1 or np.any(cp<=0):raise ValueError('invalid thermal capacity/fractions')
        measured=m*cp*(tmap-temperature0)
        sources=dict(plastic_work_J=plastic,plastic_heat_J=chi*plastic,plastic_unthermalized_J=(1-chi)*plastic,
            friction_raw_J=qpart,friction_material_heat_J=partition*qpart,friction_wall_heat_J=(1-partition)*qpart,
            contact_in_J=wall,internal_conduction_J=node)
        total=sources['plastic_heat_J']+sources['friction_material_heat_J']+wall+node;balance=measured-total
        sources.update(total_heat_input_J=total,measured_thermal_increment_J=measured,mc_delta_T_residual_J=balance)
        checks.compare('pair_heat_map_identity',thermal['pair_heat_into_i_J'],pair,1e-25,1e-12)
        checks.compare('wall_heat_boundary_once',np.asarray(boundary['contact_power_W'])*dt,wall,1e-25,1e-12)
        # Source identities use actual material journal P after its separate
        # independent identity test: final subtraction ep1-ep0 may round when
        # ep0 is nonzero. The algebraic map consumes that actual journal once.
        journal_plastic=array(mat['plastic_work_J'],(n,),'plastic journal')
        journal_friction=array(boundary['friction_raw_J'],(n,),'certified boundary raw friction')
        journal_sources=dict(sources,plastic_work_J=journal_plastic,plastic_heat_J=chi*journal_plastic,
            plastic_unthermalized_J=(1-chi)*journal_plastic,friction_raw_J=journal_friction,
            friction_material_heat_J=partition*journal_friction,friction_wall_heat_J=(1-partition)*journal_friction)
        journal_total=chi*journal_plastic+partition*journal_friction+wall+node
        journal_sources.update(total_heat_input_J=journal_total,mc_delta_T_residual_J=measured-journal_total)
        # The tiny subtraction residual is not a heat source. A relative
        # tolerance scaled to that residual is invalid under different legal
        # reduction orders; bound only this journal comparison by gamma_k
        # times the uncancelled operand magnitudes. Physical BE gates below
        # and every individual-source identity remain unchanged.
        degree=np.bincount(np.concatenate((first,second)),minlength=n)
        pair_magnitude=np.zeros(n)
        np.add.at(pair_magnitude,first,np.abs(pair));np.add.at(pair_magnitude,second,np.abs(pair))
        operation_count=64+8*friction_terms+32*degree+8*face_counts
        factor=operation_count*eps
        if np.any(factor>=1):raise ValueError('roundoff bound operation inventory unsupported')
        magnitude=(np.abs(measured)+np.abs(chi*journal_plastic)+np.abs(partition*qpart)+
            dt*(np.abs(bath)+np.abs(wallrow*tmap))+pair_magnitude)
        residual_roundoff=1e-30+factor/(1-factor)*magnitude
        global_factor=(2*n+4)*eps
        if global_factor>=1:raise ValueError('global roundoff inventory unsupported')
        global_roundoff=float(math.fsum(residual_roundoff))+global_factor/(1-global_factor)*float(np.abs(measured-journal_total).sum())
        for key,value in journal_sources.items():
            if key=='mc_delta_T_residual_J':
                checks.compare('BE_particle_'+key,thermal['particle_energy_J'][key],value,residual_roundoff)
                checks.compare('BE_ledger_'+key,thermal['ledger'][key],math.fsum(float(v) for v in value),global_roundoff)
            else:
                checks.compare('BE_particle_'+key,thermal['particle_energy_J'][key],value,1e-25,1e-12)
                checks.compare('BE_ledger_'+key,thermal['ledger'][key],math.fsum(float(v) for v in value),1e-25,1e-12)
        checks.compare('BE_wall_out',thermal['ledger']['wall_heat_out_J'],float(math.fsum((1-partition)*journal_friction)-math.fsum(wall)),1e-25,1e-12)
        # Retain the original physical BE gate's abs(total).sum() scale.
        # Uncancelled source magnitudes belong only to the roundoff diagnostic.
        scale=float(np.abs(total).sum())
        uncancelled_scale=float(np.abs(chi*plastic).sum()+np.abs(partition*qpart).sum()+np.abs(wall).sum()+np.abs(node).sum())
        checks.compare('BE_particle_balance',balance,np.zeros(n),1e-12+1e-10*scale)
        checks.compare('BE_total_balance',math.fsum(balance),0.,1e-12+1e-9*scale)
        checks.compare('BE_pair_antisymmetry',math.fsum(node),0.,1e-15+1e-10*float(np.abs(pair).sum()))
        summary.update(contact_count=len(contact_ids),contact_ids=contact_ids,feature_ids=sorted(set(features)),
            multicontact_particle_count=multi,actuator_active_particle_count=active_actuators,
            nonzero_constraint_impulse_particle_count=nonzero_reactions,contact_subintervals=subintervals,
            interval_kind_counts=kinds,events=events,negative_raw_events=negative,
            maximum_momentum_defect_Ns=maxmomentum,maximum_free_reaction_Ns=maximum_free_reaction,
            maximum_prescribed_velocity_defect_m_s=maximum_constraint_error,
            raw_friction_J=float(math.fsum(qpart)),constraint_work_J=totals['actuator_work_J'],
            contact_work_totals_J=totals,thermal_source_totals_J={k:float(math.fsum(v)) for k,v in sources.items()},
            maximum_particle_thermal_residual_J=float(np.max(abs(balance))),total_thermal_residual_J=float(math.fsum(balance)),
            thermal_absolute_total_scale_J=scale,thermal_uncancelled_source_scale_J=uncancelled_scale,
            thermal_pair_net_J=float(math.fsum(node)),
            balance_journal_roundoff_bound_J=residual_roundoff,
            balance_journal_global_roundoff_bound_J=global_roundoff,
            balance_journal_roundoff_operation_count=operation_count,
            candidate_fixed_point_convergence='NOT_TESTED_BY_MAP_AUDIT',reentry_and_parent_hash='DRIVER_RESPONSIBILITY',
            return_dissipation_is_material_heat=False,energy_residual_is_material_heat=False)
        out['passed']=all(row['passed'] for row in checks.rows.values())
        if out['passed']:out['status']='PASS_MAP_INTERNAL_ACCOUNTS_NOT_ACCEPTANCE'
    except Exception as error:
        out['error']=repr(error)
        checks.require('complete_required_structure',False,dict(error=repr(error)))
    out['checks']=list(checks.rows.values())
    summary['check_categories']=len(out['checks']);summary['checks_run']=sum(row['count'] for row in out['checks'])
    summary['failed_categories']=[row['name'] for row in out['checks'] if not row['passed']]
    return _json(out)
