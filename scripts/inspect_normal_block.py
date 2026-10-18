"""Independent saved-array audit: free-normal 3D SPH hold/lift/open benchmark.

No dynamics solver, material updater, geometry or contact implementation is
imported. Reuses only earlier independent arithmetic readers and generic VTK IO.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from inspect_rotating_deformable_block import (load,file_hash,close,array,require,
    audit_operator,reconstruct_gradient,PARAMETERS as OLD_PARAMETERS,SOURCE_PATHS as OLD_SOURCES)
from inspect_normal_contact import normal_at,tool_at,audit_normal_interval,provenance,KN,MASS,G,BULK,LAYER0,U0
from afsd_scaffold.vtk_legacy import read_point_cloud

ROOT=Path(__file__).resolve().parents[1]
CLASSIFICATION='ENGINEERING_APPROXIMATE_NORMAL_RELEASE_SPH_NOT_AFSD'
END=6e-7
SOURCES=OLD_SOURCES|{'scripts/probe_normal_deformable_block.py','src/afsd_backend/pysph/normal_penalty.py',
    'configs/benchmarks/normal_deformable_block.json','docs/NORMAL_CONTACT_CONTRACT.md',
    'third_party/PySPH/pysph/sph/rigid_body.py'}
OVERRIDES=dict(normal_force=1,layer=LAYER0,normal_stiffness=657500,penalty_factor=.05,
    contact_reference_spacing=.0005,lift_start=2e-7,lift_speed=20,end_time=6e-7,
    initial_contact_energy=16/(2*657500))


def configuration(cfg):
    require(cfg['classification']==CLASSIFICATION and cfg['production_authorized'] is False,'component policy')
    require(cfg['dt_s'] in (5e-9,2.5e-9),'registered timestep')
    require(cfg['parameters']==dict(OLD_PARAMETERS,**OVERRIDES),'literal physical parameters changed')
    require(cfg['physical_scaling']=='NONE' and cfg['tags_are_observation_only'] is True,'no scaling or active tracers')
    require(cfg['spatial_operator']=='corrected_work_conjugate','paired corrected operator')
    require(cfg['normal_law']=='UNILATERAL LINEAR PENALTY; TOP/INTERIOR XYZ FREE; BOTTOM XYZ FIXED','normal law')
    source=load(ROOT/'configs/benchmarks/normal_deformable_block.json')
    require(cfg['source_configuration']==source,'current source configuration')
    require(cfg['inherited_source_configuration']==load(ROOT/'configs/benchmarks/rotating_deformable_block.json'),'inherited physical configuration')
    require({k:v['value'] for k,v in source['parameters'].items()}==OVERRIDES,'literal source overrides')
    require(all(v.get('unit') and v.get('source') and v.get('status') for v in source['parameters'].values()),'parameter provenance')
    require(source['run_timeout_s']==170 and source['external_timeout_s']==180,'run budget')
    fields=cfg['field_status']
    require(fields['temperature']=='PRESCRIBED_INPUT_NOT_SOLVED' and fields['Es']=='DIAGNOSTIC_ONLY_NOT_CERTIFIED_FREE_ENERGY','no thermal/free-energy claim')
    require(fields['normal_force']=='SOLVED_FROM_CURRENT_OVERLAP' and fields['normal_temporal_residual']=='RAW_DIAGNOSTIC_NOT_HEAT','normal fields policy')
    return float(cfg['dt_s'])


def normal_snapshot(saved,x,ids,top,t):
    expected=normal_at(x[top],t); layer,vz=tool_at(t)
    require(saved['particle_id']==np.asarray(ids)[top].tolist(),'normal snapshot persistent IDs')
    require(saved['particle_index']==np.flatnonzero(top).tolist(),'normal snapshot index mapping')
    close(saved['time_s'],t,1e-20,'normal snapshot time')
    for key,value in [('gap_m',layer-x[top,2]),('signed_penetration_m',expected['signed']),
        ('penetration_m',expected['delta']),('normal_force_N',expected['force']),('potential_J',expected['potential'])]:
        close(saved[key],value,1e-12 if key=='normal_force_N' else 1e-17,'snapshot '+key)
    require(saved['active']==expected['active'].tolist(),'signed unilateral active state')
    close(saved['radius_m'],.00025,0.,'particle radius');close(saved['stiffness_N_m'],KN,0.,'normal stiffness')
    tool=saved['tool']; close(tool['layer'],layer,1e-17,'lifted shoulder')
    close(tool['tool_top'],.02+layer-LAYER0,1e-17,'lifted open bore top')
    close(tool['origin'],[0.,0.,layer],1e-17,'tool origin')
    close(tool['translation_velocity'],[0.,0.,vz],0.,'actual tool lift velocity')
    close(tool['omega'],[0.,0.,20.],0.,'tool omega')
    close(tool['offset'],layer-LAYER0,1e-17,'tool lift offset')
    require(np.max(expected['delta'])<=.1*.0005,'normal penetration bound')
    return expected


def history_record(r,original,ids,roles,step,dt):
    require(r['step']==step and r['particle_id']==ids and r['boundary_role']==roles,'history identity/roles')
    require(r['tracer_type']==[pid%3 for pid in ids] and r['as_rs_tag']==[0]*64,'observational labels')
    t=step*dt;close(r['time_s'],t,1e-20,'history time')
    x=array(r['position_m'],(64,3),'position');v=array(r['velocity_m_s'],(64,3),'velocity')
    mass=array(r['mass_kg'],(64,),'mass');rho=array(r['density_kg_m3'],(64,),'density')
    close(mass,np.full(64,MASS),1e-21,'unscaled masses');close(mass.sum(),2.16e-5,2.16e-17,'total mass')
    require(np.max(np.abs(rho/2700-1))<=.01,'density bound')
    require(np.max(np.linalg.norm(v,axis=1))<=1.,'material speed bound')
    require(np.max(np.linalg.norm(x-original,axis=1))<=1e-5,'material displacement bound')
    fixed=np.asarray(roles)==1;top=original[:,2]==original[:,2].max()
    close(x[fixed],original[fixed],2e-17,'bottom fixed position');close(v[fixed],np.zeros((16,3)),1e-12,'bottom fixed velocity')
    s=array(r['deviatoric_stress_Pa'],(64,3,3),'deviatoric stress')
    close(s,s.swapaxes(1,2),1e-6,'symmetric stress');close(np.trace(s,axis1=1,axis2=2),np.zeros(64),1e-6,'deviatoric trace')
    pressure=BULK*(rho/2700-1);close(r['pressure_Pa'],pressure,1e-6,'EOS pressure')
    close(r['stress_Pa'],s-pressure[:,None,None]*np.eye(3),1e-6,'total stress')
    close(r['equivalent_stress_Pa'],np.sqrt(1.5*np.sum(s*s,axis=(1,2))),1e-6,'equivalent stress')
    ep=array(r['ep'],(64,),'plastic strain');require(np.min(ep)>=0,'plastic strain sign')
    normal=normal_snapshot(r['normal_contact'],x,ids,top,t);potential=float(normal['potential'].sum())
    require(r['contact_count']==int(normal['active'].sum()),'actual active contact count')
    close(r['normal_potential_J'],potential,1e-17,'actual contact potential')
    close(r['initial_normal_potential_J'],U0,1e-17,'supplied initial contact energy')
    close(r['max_penetration_m'],normal['delta'].max(),1e-17,'maximum penetration')
    layer,vz=tool_at(t);close(r['tool_origin_m'],[0.,0.,layer],1e-17,'record tool origin')
    close(r['tool_translation_velocity_m_s'],[0.,0.,vz],0.,'record tool velocity')
    close(r['tool_angle_rad'],20*t,1e-20,'record tool angle')
    kinetic=float(.5*np.sum(mass[:,None]*v*v))
    shear=float(np.sum(mass/rho*np.sum(s*s,axis=(1,2))/(4*G)))
    ratio=rho/2700;volume=float(np.sum(mass*BULK/2700*(np.log(ratio)+1/ratio-1)))
    close(r['kinetic_J'],kinetic,1e-18,'actual kinetic energy')
    close(r['Es_diagnostic_J'],shear+volume,1e-16,'diagnostic Es not free energy')
    c=r['cumulative'];residual=kinetic+shear+volume+c['plastic_work_J']+c['friction_raw_J']+potential-U0-c['wall_work_J']-c['constraint_work_J']
    close(r['energy_diagnostic_residual_J'],residual,1e-16,'reported global energy diagnostic')
    warning=abs(residual)>1e-12+.01*(U0+abs(c['wall_work_J'])+abs(c['constraint_work_J'])+abs(c['friction_raw_J']))
    require(r['energy_status']==('ENERGY_WARNING' if warning else 'DIAGNOSTIC_WITHIN_SCREEN_NOT_CERTIFIED_FREE_ENERGY'),'energy warning label')
    return dict(x=x,v=v,mass=mass,rho=rho,s=s,ep=ep,K=kinetic,U=potential,normal=normal,residual=residual,warning=warning)


def contact_trace(trace,index,pid,role,x,vb,va,force,dt,step):
    require(trace['particle_index']==index and trace['particle_id']==pid and trace['step']==step,'kick trace identity')
    require(trace['stage']=='MIDPOINT_FROZEN_GEOMETRY_KICK','kick stage')
    t=(step-.5)*dt;close(trace['time_mid_s'],t,1e-20,'contact midpoint time')
    close(trace['position_m'],x,2e-17,'contact application position')
    close(trace['material_velocity_before_m_s'],vb,1e-12,'trace before velocity')
    close(trace['internal_force_N'],force,1e-12,'trace internal force');close(trace['mass_kg'],MASS,1e-21,'trace mass')
    free=[0.,0.,0.] if role==1 else [1.,1.,1.];require(trace['free_axes']==free,'top/interior free mobility')
    normal=normal_at(np.array([x]),t);active=bool(normal['active'][0]);fn=float(normal['force'][0])
    faces=trace['contact_rows'];require(len(faces)==int(active),'active contact inventory')
    layer,vz=tool_at(t);wall=np.array([[-20*x[1],20*x[0],vz]]) if active else np.empty((0,3))
    for face in faces:
        require(face['particle_id']==pid and face['feature_id']==1 and face['contact_id']==str(pid)+':1','active shoulder feature')
        close(face['application_point'],x,2e-17,'force application')
        close(face['contact_point'],[x[0],x[1],layer],2e-17,'face projection')
        close(face['normal'],[0.,0.,-1.],0.,'normal direction')
        close(face['gap'],layer-x[2],2e-17,'face gap');close(face['normal_force_N'],fn,1e-12,'actual gap-derived Fn')
        close(face['penetration'],normal['delta'][0],2e-17,'active overlap')
        close(face['wall_velocity'],wall[0],1e-12,'actual moving wall velocity')
        close(face['wall_velocity_at_face'],wall[0],1e-12,'shoulder diagnostic face velocity')
        close(face['translation_velocity'],[0.,0.,vz],0.,'lift translation at contact')
        close(face['body_origin'],[0.,0.,layer],2e-17,'moving contact origin')
        close(face['omega'],[0.,0.,20.],0.,'contact angular velocity')
    result=trace['result'];close(result['constraint_projector'],np.diag(free),0.,'free constraint projector')
    close(result['velocity'],va,1e-12,'kick endpoint')
    require(bool(result['intervals']),'kick accepted intervals')
    values=dict(internal=0.,normal=0.,normal_wall=0.,wall=0.,constraint=0.,kinetic=0.,q=0.,negative=0)
    external=np.zeros(3);contact=np.zeros(3);workv=np.zeros(3);elapsed=0.;previous=vb;minq=None
    for sub in result['intervals']:
        a=audit_normal_interval(sub,previous,force,wall,fn,free,MASS,dt,elapsed,vb)
        previous=a['v_after'];elapsed+=a['h'];workv+=a['work_velocity']
        external+=a['external_impulse'];contact+=a['contact_impulse']
        for key in values:values[key]+=a[key]
        if a['min_q'] is not None:minq=a['min_q'] if minq is None else min(minq,a['min_q'])
    close(elapsed,dt,1e-20,'whole kick duration');close(previous,va,1e-12,'whole kick final velocity')
    for target,source in [('other_work_J','internal'),('wall_work_J','wall'),('constraint_work_J','constraint'),('delta_kinetic_J','kinetic')]:
        close(result[target],values[source],1e-15,'kick summed '+target)
    close(np.sum(result['friction_heat_J']),values['q'],1e-18,'summed raw friction')
    close(result['kinetic_work_residual_J'],0.,1e-15,'kinetic work diagnostic')
    return dict(**values,external=external,contact=contact,workv=workv/dt,min_q=minq)


def vtk_record(path,r,original):
    v=read_point_cloud(path);scalars=v['scalars'];step=r['step']
    require(CLASSIFICATION in v['title'] and 'NO_THERMAL_SOLUTION' in v['title'],'VTK honest classification')
    names={'particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','boundary_role','as_rs_tag',
        'mass_kg','density_kg_per_m3','pressure_Pa','equivalent_stress_Pa','equivalent_plastic_strain','smoothing_length_m',
        'signed_gap_m','normal_signed_penetration_m','normal_penetration_m','normal_force_N','normal_potential_J','contact_active',
        *('stress_'+name+'_Pa' for name in ('xx','yy','zz','xy','xz','yz'))}
    if step:names|={'equivalent_deviatoric_strain_rate_per_s','equivalent_plastic_strain_rate_per_s'}
    require(set(scalars)==names and set(v['vectors'])=={'velocity_m_per_s','initial_position_m'},'VTK omitted/fabricated fields')
    close(v['points'],r['position_m'],0.,'VTK positions');close(v['vectors']['velocity_m_per_s'],r['velocity_m_s'],0.,'VTK velocity')
    close(v['vectors']['initial_position_m'],original,0.,'VTK initial geometry')
    for key in ('particle_id','tracer_type','boundary_role','as_rs_tag'):require(scalars[key]==r[key],'VTK '+key)
    fields=dict(mass_kg='mass_kg',density_kg_per_m3='density_kg_m3',pressure_Pa='pressure_Pa',equivalent_stress_Pa='equivalent_stress_Pa',equivalent_plastic_strain='ep')
    if step:fields.update(equivalent_deviatoric_strain_rate_per_s='strain_rate_per_s',equivalent_plastic_strain_rate_per_s='ep_rate_per_s')
    for target,source in fields.items():close(scalars[target],r[source],0.,'VTK '+target)
    for a,b,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        close(scalars['stress_'+name+'_Pa'],np.asarray(r['stress_Pa'])[:,a,b],0.,'VTK '+name)
    normal=normal_at(r['position_m'],r['time_s']);x=np.asarray(r['position_m'])
    for key,value in [('signed_gap_m',normal['layer']-x[:,2]),('normal_signed_penetration_m',normal['signed']),
        ('normal_penetration_m',normal['delta']),('normal_force_N',normal['force']),('normal_potential_J',normal['potential']),('contact_active',normal['active'].astype(int))]:
        close(scalars[key],value,1e-12 if key=='normal_force_N' else 1e-17,'VTK '+key)
    close(scalars['material_id'],np.ones(64),0.,'unchanged tracer material')
    close(scalars['smoothing_length_m'],np.full(64,.00065),1e-18,'unscaled support')
    for key in ('component_id','wholebar_tag','actuator'):close(scalars[key],np.zeros(64),0.,'VTK '+key)


def audit_block(directory):
    directory=Path(directory);m=provenance(directory,'COMPLETE_NORMAL_DEFORMABLE_COMPONENT_NOT_AFSD')
    dt=configuration(m['configuration']);count=round(END/dt)
    require(m['mode']=='FULL_REGISTERED_NORMAL_COMPONENT' and set(m['source_sha256'])==SOURCES,'scope/source inventory')
    expected={'history.json','steps.jsonl','last_accepted.npz','prelift_accepted.npz',*(f'frame_{k:04d}.vtk' for k in range(count+1))}
    require(set(m['artifact_sha256'])==expected,'exact complete artifact inventory')
    require(m['accepted_steps']==m['target_steps']==count,'complete registered step count')
    close(m['accepted_time_s'],END,1e-20,'common physical endpoint')
    require(m['maximum_wall_seconds']<=170 and m['wall_seconds']<=170,'bounded wall time')
    close(m['initial_normal_potential_J'],U0,1e-17,'initial potential source')
    h=load(directory/'history.json');journal=[json.loads(line) for line in (directory/'steps.jsonl').read_text(encoding='utf-8').splitlines()]
    require(len(h)==count+1 and len(journal)==count,'history/journal completeness')
    original=np.array(np.meshgrid(.01+(np.arange(4)-1.5)*.0005,(np.arange(4)-1.5)*.0005,(np.arange(4)+.5)*.0005,indexing='ij')).reshape(3,-1).T
    ids=list(range(1,65));top=original[:,2]==.00175;bottom=original[:,2]==.00025
    roles=np.where(bottom,1,np.where(top,2,0)).tolist();state=history_record(h[0],original,ids,roles,0,dt)
    for actual,expected,name in [(state['x'],original,'geometry'),(state['v'],np.zeros((64,3)),'velocity'),
        (state['s'],np.zeros((64,3,3)),'stress'),(state['rho'],np.full(64,2700.),'density'),(state['ep'],np.zeros(64),'plastic strain')]:
        close(actual,expected,0.,'initial '+name)
    vtk_record(directory/'frame_0000.vtk',h[0],original)
    sums={key:0. for key in h[0]['cumulative']};require(all(value==0 for value in h[0]['cumulative'].values()),'initial ledger zero')
    maxcondition=0.;maxpower=0.;maxforce=0.;maxenergy=0.;warnings=0;negative=0;minq=None
    external_total=np.zeros(3);angular_total=np.zeros(3);opening={};maximum_inner_vz=0.
    for step,(row,r) in enumerate(zip(journal,h[1:]),1):
        require(row['step']==step and row['particle_id']==ids and row['boundary_role']==roles,'journal identity')
        require(row['stage']=='ACCEPTED_DRIFT_KICK_DRIFT_COMPLETE','only accepted state journals')
        close(row['dt_s'],dt,0.,'journal dt');close(row['time_mid_s'],(step-.5)*dt,1e-20,'journal midpoint time')
        current=history_record(r,original,ids,roles,step,dt)
        xm=array(row['x_mid_m'],(64,3),'drift midpoint');vb=array(row['v_before_m_s'],(64,3),'pre-kick velocity');va=array(row['v_after_m_s'],(64,3),'post-kick velocity')
        close(vb,state['v'],1e-12,'history kick input');close(va,current['v'],1e-12,'history kick output')
        close(xm,state['x']+.5*dt*vb,2e-17,'first drift');close(current['x'],xm+.5*dt*va,2e-17,'second drift')
        close(row['x_before_m'],state['x'],0.,'saved pre-step x');close(row['x_after_m'],current['x'],0.,'saved post-step x')
        rho=array(row['rho_before_kg_m3'],(64,),'operator density');close(rho,state['rho'],1e-12,'density chronology')
        close(row['force_density_kg_m3'],rho,0.,'same force and moment volume');close(row['mass_kg'],current['mass'],0.,'journal mass')
        i,j=np.asarray(row['i'],int),np.asarray(row['j'],int);g=np.asarray(row['g']);B=np.asarray(row['B'])
        sigma=array(row['sigma_mid_Pa'],(64,3,3),'force stress');force=array(row['internal_force_N'],(64,3),'internal force')
        L=reconstruct_gradient(vb,i,j,g,current['mass']/rho,B);close(row['L_before_per_s'],L,1e-10,'precontact gradient')
        drho=-rho*np.trace(L,axis1=1,axis2=2);close(row['drho_kg_m3_s'],drho,1e-6,'continuity rate')
        close(current['rho'],rho+dt*drho,1e-10,'density integration');close(row['rho_after_kg_m3'],current['rho'],0.,'final density')
        close(row['rho_mid_kg_m3'],rho+.5*dt*drho,1e-10,'midpoint EOS density')
        close(row['S_before_Pa'],state['s'],1e-6,'old material stress');close(row['S_after_Pa'],current['s'],1e-6,'accepted material stress')
        close(sigma,.5*(state['s']+current['s'])-BULK*((rho+.5*dt*drho)/2700-1)[:,None,None]*np.eye(3),1e-6,'stress used by internal force')
        acoustic=.1*.00065/(np.sqrt((BULK+4*G/3)/2700)+np.max(np.linalg.norm(vb,axis=1)))
        require(dt<=min(acoustic,.1*np.sqrt(MASS/KN)),'acoustic and contact timestep bounds')
        for name,x,time in [('normal_before',state['x'],(step-1)*dt),('normal_midpoint',xm,(step-.5)*dt),('normal_after',current['x'],step*dt)]:
            normal_snapshot(row[name],x,ids,top,time)
        traces=row['traces'];require(len(traces)==64,'all material kick traces')
        workv=np.zeros((64,3));external=np.zeros((64,3));contact=np.zeros((64,3));inc=dict.fromkeys(sums,0.)
        for index,trace in enumerate(traces):
            a=contact_trace(trace,index,ids[index],roles[index],xm[index],vb[index],va[index],force[index],dt,step)
            workv[index]=a['workv'];external[index]=a['external'];contact[index]=a['contact']
            for target,source in [('friction_raw_J','q'),('wall_work_J','wall'),('constraint_work_J','constraint'),
                ('internal_work_J','internal'),('normal_work_J','normal'),('normal_wall_work_J','normal_wall'),('delta_kinetic_J','kinetic')]:inc[target]+=a[source]
            negative+=a['negative']
            if a['min_q'] is not None:minq=a['min_q'] if minq is None else min(minq,a['min_q'])
        close(row['work_velocity_m_s'],workv,1e-12,'actual subinterval work velocity')
        close(row['endpoint_midpoint_velocity_m_s'],.5*(vb+va),1e-12,'endpoint midpoint separately recorded')
        op=audit_operator(xm,workv,current['mass'],rho,i,j,g,B,np.asarray(row['L_midpoint_per_s']),force,sigma)
        maxcondition=max(maxcondition,op['maximum_condition']);maxpower=max(maxpower,op['power_residual_W'])
        maxforce=max(maxforce,float(np.linalg.norm(force,axis=1).max()))
        inc['material_old_work_J']=float(dt*np.einsum('n,nij,nij',current['mass']/rho,sigma,L))
        inc['midpoint_stress_work_J']=dt*op['stress_power_W']
        inc['plastic_work_J']=float(np.sum(row['plastic_work_J']));inc['return_dissipation_J']=float(np.sum(row['return_dissipation_J']))
        # This short gate is expected elastic; if it leaves that domain, do not
        # grant unbenchmarked plastic flow validation from a readback identity.
        close(current['ep'],np.zeros(64),0.,'registered short elastic domain')
        close(row['plastic_work_J'],np.zeros(64),0.,'actual zero plastic work')
        close(row['return_dissipation_J'],np.zeros(64),0.,'actual zero return dissipation')
        D=.5*(L+L.swapaxes(1,2));D-=np.trace(L,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
        close(r['strain_rate_per_s'],np.sqrt(2/3*np.sum(D*D,axis=(1,2))),1e-10,'material strain rate time layer')
        close(r['ep_rate_per_s'],np.zeros(64),0.,'elastic ep rate')
        inc['normal_relative_work_J']=inc['normal_work_J']-inc['normal_wall_work_J']
        du=current['U']-state['U'];inc['normal_temporal_residual_J']=du+inc['normal_relative_work_J']
        for name,value in [('potential_change_J',du),('material_work_J',inc['normal_work_J']),('wall_work_J',inc['normal_wall_work_J']),('temporal_residual_J',inc['normal_temporal_residual_J'])]:
            close(row['normal_energy'][name],value,1e-15,'separate normal ledger '+name)
        require(row['normal_energy']['status']=='RAW_TEMPORAL_DISCRETIZATION_DIAGNOSTIC_NOT_HEAT','normal residual not heat')
        close(inc['internal_work_J']+inc['midpoint_stress_work_J'],0.,1e-15,'actual spatial work identity')
        dp=np.sum(current['mass'][:,None]*(va-vb),axis=0)
        dl=np.sum(current['mass'][:,None]*(np.cross(current['x'],va)-np.cross(state['x'],vb)),axis=0)
        angular=np.cross(xm,external).sum(axis=0)
        close(dp,external.sum(axis=0),1e-17,'whole body external impulse');close(dl,angular,1e-19,'whole body angular impulse')
        wall_vz=tool_at((step-.5)*dt)[1]
        close(inc['wall_work_J'],20*np.cross(xm,contact).sum(axis=0)[2]+wall_vz*contact[:,2].sum(),1e-15,'rotation plus lift conjugate work')
        external_total+=external.sum(axis=0);angular_total+=angular
        for key,value in inc.items():
            close(row['increments'][key],value,1e-15,'saved increment '+key)
            sums[key]+=value;close(r['cumulative'][key],sums[key],1e-14,'cumulative '+key)
        kinetic_residual=current['K']+sums['friction_raw_J']-sums['wall_work_J']-sums['constraint_work_J']-sums['internal_work_J']-sums['normal_relative_work_J']
        close(kinetic_residual,0.,1e-15,'global contact normal kinetic balance');close(r['kinetic_contact_residual_J'],kinetic_residual,1e-15,'saved kinetic residual')
        if step==1:
            wall=np.column_stack((-20*original[top,1],20*original[top,0],np.zeros(16)))
            impulse=.2*dt*wall/np.linalg.norm(wall,axis=1)[:,None];impulse[:,2]=-dt
            expected_v=np.zeros((64,3));expected_v[top]=impulse/MASS
            close(va,expected_v,1e-12,'independent first free-normal velocity')
            close(current['x'],original+.5*dt*expected_v,2e-17,'independent first position')
            q=float(np.sum(impulse[:,:2]*wall[:,:2])-np.sum(impulse[:,:2]**2)/(2*MASS))
            close(sums['friction_raw_J'],q,1e-18,'independent first signed friction')
        for index,pid in enumerate(np.asarray(ids)[top]):
            if state['normal']['active'][index] and not current['normal']['active'][index] and str(pid) not in opening:
                opening[str(pid)]=[(step-1)*dt,step*dt]
        maxenergy=max(maxenergy,abs(current['residual']));warnings+=int(current['warning'])
        maximum_inner_vz=max(maximum_inner_vz,float(np.max(np.abs(va[np.asarray(roles)==0,2]))))
        vtk_record(directory/f'frame_{step:04d}.vtk',r,original);state=current
    require(len(opening)==16 and not np.any(state['normal']['active']),'all top particles actually opened')
    close(state['normal']['force'],np.zeros(16),0.,'open force exactly zero');close(state['U'],0.,0.,'open potential exactly zero')
    require(maxforce>1e-6 and maximum_inner_vz>1e-8,'nonzero free interior dynamics')
    require(np.max(np.abs(state['s'][:,0,2]))>1 and np.max(np.abs(state['s'][:,1,2]))>1,'nontrivial xz and yz stress')
    require(m['final_record']==h[-1],'manifest accepted history consistency')
    return dict(status='PASS_NORMAL_BLOCK_READBACK_NOT_AFSD',classification=CLASSIFICATION,production_authorized=False,
        case=str(directory),dt_s=dt,steps=count,end_time_s=END,frames_checked=count+1,source_sha256=m['source_sha256'],
        manifest_sha256=file_hash(directory/'manifest.json'),reader_sha256=file_hash(Path(__file__)),
        reader_dependencies_sha256={p:file_hash(ROOT/'scripts'/p) for p in ('inspect_normal_contact.py','inspect_rotating_deformable_block.py')},
        minimum_neighbor_rank=3,maximum_condition=maxcondition,maximum_power_residual_W=maxpower,maximum_internal_force_N=maxforce,
        maximum_inner_z_velocity_m_s=maximum_inner_vz,cumulative=sums,minimum_raw_friction_J=minq,negative_raw_intervals=negative,
        initial_contact_energy_J=U0,final_contact_energy_J=state['U'],final_active_contacts=0,opening_brackets_by_id_s=opening,
        maximum_energy_diagnostic_residual_J=maxenergy,energy_warning_frames=warnings,energy_certified=False,
        final_energy_diagnostic_residual_J=state['residual'],final_ep_max=float(state['ep'].max()),
        total_external_impulse_Ns=external_total.tolist(),total_external_angular_impulse_Nms=angular_total.tolist(),
        checkpoints='HASH_CHECKED; accepted/prelift loaders and continuation checked separately, not certified by file hash',
        continuum_reference='NOT_AVAILABLE; independent first-step and algebra only, no exact whole-trajectory solution')


def compare_blocks(coarse,fine):
    a,b=audit_block(coarse),audit_block(fine)
    require(a['dt_s']==5e-9 and b['dt_s']==2.5e-9 and a['end_time_s']==b['end_time_s']==END,'same registered common endpoint')
    require(a['source_sha256']==b['source_sha256'],'same physical source')
    ca=load(Path(coarse)/'manifest.json')['configuration'];cb=load(Path(fine)/'manifest.json')['configuration']
    require({k:v for k,v in ca.items() if k!='dt_s'}=={k:v for k,v in cb.items() if k!='dt_s'},'same physical parameters geometry boundaries')
    ha,hb=load(Path(coarse)/'history.json'),load(Path(fine)/'history.json');metrics={}
    for field,limit in dict(position_m=1e-7,velocity_m_s=1e-2,stress_Pa=1e5,density_kg_m3=1.).items():
        differences=[]
        for k,r in enumerate(ha):
            require(r['particle_id']==hb[2*k]['particle_id'],'common time persistent IDs')
            close(r['time_s'],hb[2*k]['time_s'],1e-20,'common physical frame')
            differences.append(float(np.max(np.abs(np.asarray(r[field])-np.asarray(hb[2*k][field])))))
        metrics[field]=dict(maximum_common_history_difference=max(differences),endpoint_difference=differences[-1],limit=limit,within_screen=max(differences)<=limit)
    for field in ('friction_raw_J','wall_work_J'):
        differences=[abs(r['cumulative'][field]-hb[2*k]['cumulative'][field]) for k,r in enumerate(ha)]
        metrics[field]=dict(maximum_common_history_difference=max(differences),endpoint_difference=differences[-1],limit=1e-7,within_screen=max(differences)<=1e-7)
    opening={}
    for pid,bracket in a['opening_brackets_by_id_s'].items():
        finebracket=b['opening_brackets_by_id_s'][pid];difference=abs(bracket[1]-finebracket[1])
        opening[pid]=dict(coarse=bracket,fine=finebracket,upper_sample_difference_s=difference,within_screen=difference<=1e-8)
    within=all(v['within_screen'] for v in metrics.values()) and all(v['within_screen'] for v in opening.values())
    return dict(status='NORMAL_BLOCK_COMPARISON_COMPLETE_NOT_CONVERGENCE_CERTIFICATE',production_authorized=False,
        screening_status='WITHIN_SCREEN' if within else 'SENSITIVE_WARNING',common_end_s=END,common_frames=len(ha),
        metrics=metrics,opening_comparison_by_id=opening,cases=[a,b],
        note='Fixed physical problem, 121 common samples; no measured convergence order or paper-scale validation.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block',type=Path);parser.add_argument('--coarse',type=Path);parser.add_argument('--fine',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    require(not args.output.exists(),'fresh reader output required')
    if bool(args.block)==bool(args.coarse or args.fine) or (not args.block and not(args.coarse and args.fine)):
        parser.error('use --block OR --coarse and --fine')
    try:
        result=audit_block(args.block) if args.block else compare_blocks(args.coarse,args.fine);code=0
    except Exception as error:
        result=dict(status='FAIL_NORMAL_BLOCK_READBACK',error=str(error),production_authorized=False,reader_sha256=file_hash(Path(__file__)));code=2
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,allow_nan=False));return code


if __name__=='__main__':raise SystemExit(main())
