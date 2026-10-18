"""Independent array audit for initial-material axial drive and free release.

Only independent audit arithmetic and generic VTK IO are reused. No solver,
actuator implementation, material updater, or contact reference generator import.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from inspect_rotating_deformable_block import (load,file_hash,close,array,require,canonical_hash,
    audit_operator,reconstruct_gradient,PARAMETERS as OLD_PARAMETERS,SOURCE_PATHS as OLD_SOURCES)
from afsd_scaffold.vtk_legacy import read_point_cloud

ROOT=Path(__file__).resolve().parents[1]
CLASSIFICATION='ENGINEERING_APPROXIMATE_LOCAL_AXIAL_FEED_NOT_AFSD'
MASS=2700*.0005**3;KN=657500.;G=26.3e9;BULK=2*G*1.33/(3*.34)
END=8e-7;RELEASE=4e-7;SPEED=.00212;K0=1.213488e-11;P0=np.array([0.,0.,-1.1448e-8])
OVERRIDES=dict(omega=0,traverse=0,friction=0,layer=-.009,feed=.00212,actuator_depth=.0005,
    marker_low=-.01,marker_high=-.008,normal_stiffness=657500,release_time=4e-7,end_time=8e-7)


def literal_geometry():
    return np.array(np.meshgrid((np.arange(4)-1.5)*.0005,(np.arange(4)-1.5)*.0005,
        -.01+(np.arange(4)+.5)*.0005,indexing='ij')).reshape(3,-1).T


def actuator_reference(z0,time):
    """Literal fixed material subset; current position is intentionally absent."""
    eligible=(np.asarray(z0)>-.0085)&(np.asarray(z0)<=-.008)
    active=eligible&(time<RELEASE)
    free=np.ones((len(z0),3));free[active,2]=0.
    prescribed=np.zeros_like(free);prescribed[active,2]=-SPEED
    return dict(eligible=eligible,active=active,free_axes=free,prescribed_velocity=prescribed,
        phase='DRIVEN' if time<RELEASE else 'RELEASED',commanded_displacement_m=-SPEED*min(time,RELEASE))


def normal_reference(x):
    x=np.asarray(x);gap=x[:,2]+.01;signed=.00025-gap;active=signed>0
    overlap=np.where(active,signed,0.)
    return dict(gap=gap,signed=signed,overlap=overlap,active=active,force=KN*overlap,potential=.5*KN*overlap**2)


def audit_control(saved,z0,t):
    expected=actuator_reference(z0,t)
    for key in ('eligible','active'):
        require(np.asarray(saved[key]).tolist()==expected[key].tolist(),'actuator original-material '+key)
    for key in ('free_axes','prescribed_velocity'):
        close(saved[key],expected[key],0.,'actuator '+key)
    require(saved['phase']==expected['phase'],'actuator phase')
    close(saved['commanded_displacement_m'],expected['commanded_displacement_m'],2e-17,'actuator command not free path')
    require(saved['displacement_status']=='ACTUATOR_COMMAND_ONLY_NOT_FREE_MATERIAL_PATH','command provenance')
    return expected


def audit_interval(a,vbefore,force,fn,free,prescribed,mass,dt,start):
    """Independent +z base, mu=0 impulse/work audit; no Q sign repair."""
    h=float(a['dt_s']);require(np.isfinite(h) and 0<h<=dt,'subinterval duration')
    close(a['start_s'],start,1e-20,'subinterval chronology')
    vb=array(a['v_before'],(3,),'before velocity');va=array(a['v_after'],(3,),'after velocity')
    close(vb,vbefore,1e-12,'subinterval velocity continuity')
    nc=1 if fn>0 else 0
    jn=np.asarray(a['normal_impulse']).reshape(nc,3);jt=np.asarray(a['tangential_impulse']).reshape(nc,3)
    reaction=array(a['constraint_impulse'],(3,),'actuator impulse');free=np.asarray(free);prescribed=np.asarray(prescribed)
    close(jn,np.tile([0.,0.,fn*h],(nc,1)),1e-18,'gap-derived base impulse')
    close(jt,np.zeros((nc,3)),0.,'explicit zero-mu impulse')
    close(a['friction_heat_raw_J'],np.zeros(nc),0.,'explicit zero-mu raw heat')
    close(a['friction_rounding_bound_J'],np.zeros(nc),0.,'zero impulse rounding bound')
    close(mass*(va-vb),h*np.asarray(force)+jn.sum(axis=0)+reaction,1e-17,'actual impulse balance')
    close(free*reaction,np.zeros(3),1e-17,'no reaction on free DOF')
    close((1-free)*(va-prescribed),np.zeros(3),1e-12,'prescribed axial velocity only')
    mid=.5*(vb+va);tangent=np.tile(mid,(nc,1));tangent[:,2]=0.
    close(np.asarray(a['tangential_relative_midpoint']).reshape(nc,3),tangent,1e-12,'actual tangential midpoint')
    internal=float(h*np.asarray(force)@mid);normal=float(jn.sum(axis=0)@mid);actuator=float(reaction@mid)
    kinetic=.5*mass*float(va@va-vb@vb);residual=kinetic-internal-normal-actuator
    tol=1e-18+1e-9*(abs(kinetic)+abs(internal)+abs(normal)+abs(actuator))
    close(residual,0.,tol,'normal-actuator kinetic work identity')
    for key,value in [('other_work_J',internal),('constraint_work_J',actuator),('wall_work_J',0.),('delta_kinetic_J',kinetic)]:
        close(a[key],value,tol,'saved '+key)
    close(np.sum(a['normal_work_J']),normal,tol,'normal material work')
    close(a['kinetic_work_residual_J'],residual,tol,'saved local kinetic residual')
    return dict(v_after=va,h=h,work_velocity=h*mid,internal=internal,normal=normal,actuator=actuator,
        kinetic=kinetic,normal_impulse=jn.sum(axis=0),actuator_impulse=reaction,external=jn.sum(axis=0)+reaction)


def configuration(cfg):
    require(cfg['classification']==CLASSIFICATION and cfg['production_authorized'] is False,'configuration scope')
    expected=dict(OLD_PARAMETERS,**OVERRIDES);expected.pop('normal_force')
    require(cfg['parameters']==expected,'frozen physical parameter literals')
    require(cfg['dt_s'] in (5e-9,2.5e-9),'registered timestep')
    require(cfg['physical_scaling']=='NONE' and cfg['tags_are_observation_only'] is True,'no scaling or active tracer')
    require(cfg['spatial_operator']=='corrected_work_conjugate','paired corrected operator')
    source=load(ROOT/'configs/benchmarks/axial_feed.json')
    require(cfg['source_configuration']==source,'source configuration matches disk')
    require(cfg['inherited_source_configuration']==load(ROOT/'configs/benchmarks/rotating_deformable_block.json'),'inherited source configuration')
    require({k:v['value'] for k,v in source['parameters'].items()}==OVERRIDES,'source literal values')
    require(source['remove_inherited_parameters']==['normal_force'],'no inherited prescribed normal load')
    require(all(v.get('unit') and v.get('status') and v.get('source') for v in source['parameters'].values()),'parameter units/provenance')
    require(source['run_timeout_s']==170 and source['external_timeout_s']==180,'registered wall limit')
    status=cfg['field_status'];require(status['temperature']=='PRESCRIBED_INPUT_NOT_SOLVED' and
        status['Es']=='DIAGNOSTIC_ONLY_NOT_CERTIFIED_FREE_ENERGY','field certification limits')
    return float(cfg['dt_s'])


def audit_normal_snapshot(saved,x,ids,t):
    expected=normal_reference(x)
    require(saved['particle_id']==ids and saved['particle_index']==list(range(64)),'base snapshot all-particle identity')
    close(saved['time_s'],t,1e-20,'base snapshot time')
    for key,value in [('gap_m',expected['gap']),('signed_penetration_m',expected['signed']),
        ('penetration_m',expected['overlap']),('normal_force_N',expected['force']),('potential_J',expected['potential'])]:
        close(saved[key],value,1e-12 if key=='normal_force_N' else 1e-17,'base snapshot '+key)
    require(saved['active']==expected['active'].tolist(),'actual unilateral base state')
    require(np.max(expected['overlap'])<=.00005,'normal penetration limit')
    return expected


def history_record(r,original,ids,step,dt):
    t=step*dt
    require(r['step']==step and r['particle_id']==ids,'record persistent identity')
    expected_tags=(np.max(np.abs(original[:,:2]),axis=1)>=.001-.0005*.51).astype(int).tolist()
    require(r['tracer_type']==expected_tags and r['as_rs_tag']==[0]*64,'initial observational labels')
    require(r['boundary_role']==(2*(original[:,2]>-.0085)).tolist(),'initial-material role metadata')
    close(r['time_s'],t,1e-20,'record time')
    x=array(r['position_m'],(64,3),'position');v=array(r['velocity_m_s'],(64,3),'velocity')
    m=array(r['mass_kg'],(64,),'mass');rho=array(r['density_kg_m3'],(64,),'density')
    close(m,np.full(64,MASS),1e-21,'unscaled physical masses');close(m.sum(),2.16e-5,2.16e-17,'mass conservation')
    require(np.max(np.abs(rho/2700-1))<=.01,'density limit')
    require(np.max(np.linalg.norm(v,axis=1))<=1.,'speed limit')
    require(np.max(np.linalg.norm(x-original,axis=1))<=1e-5,'displacement limit')
    control=audit_control(r['control'],original[:,2],t)
    if t<=RELEASE:
        selected=control['eligible']
        close(v[selected,2],np.full(16,-SPEED),1e-12,'driven original-layer z speed')
        close(x[selected,2],original[selected,2]-SPEED*t,2e-17,'driven axial displacement')
    s=array(r['deviatoric_stress_Pa'],(64,3,3),'deviatoric stress')
    close(s,s.swapaxes(1,2),1e-6,'symmetric stress');close(np.trace(s,axis1=1,axis2=2),np.zeros(64),1e-6,'deviatoric trace')
    pressure=BULK*(rho/2700-1);close(r['pressure_Pa'],pressure,1e-6,'EOS pressure')
    close(r['stress_Pa'],s-pressure[:,None,None]*np.eye(3),1e-6,'total stress')
    close(r['equivalent_stress_Pa'],np.sqrt(1.5*np.sum(s*s,axis=(1,2))),1e-6,'equivalent stress')
    ep=array(r['ep'],(64,),'plastic strain');require(np.min(ep)>=0,'nonnegative plastic strain')
    normal=audit_normal_snapshot(r['normal_contact'],x,ids,t);potential=float(normal['potential'].sum())
    initial_u=float(normal_reference(original)['potential'].sum())
    close(r['normal_potential_J'],potential,1e-18,'base potential')
    close(r['initial_normal_potential_J'],initial_u,1e-25,'initial potential is not fabricated load')
    require(r['contact_count']==int(normal['active'].sum()),'actual base contact count')
    close(r['max_penetration_m'],normal['overlap'].max(),1e-17,'reported maximum penetration')
    close(r['initial_kinetic_J'],K0,1e-25,'explicit initial kinetic energy')
    close(r['initial_linear_momentum_kg_m_s'],P0,1e-22,'explicit initial momentum')
    v0=np.zeros((64,3));v0[original[:,2]>-.0085,2]=-SPEED
    initial_L=np.sum(MASS*np.cross(original,v0),axis=0)
    close(r['initial_angular_momentum_kg_m2_s'],initial_L,1e-24,'explicit initial angular momentum')
    close(r['tool_origin_m'],[0.,0.,-.009],0.,'fixed open tool origin')
    close(r['tool_angle_rad'],0.,0.,'explicit nonrotating probe')
    kinetic=float(.5*np.sum(m[:,None]*v*v));shear=float(np.sum(m/rho*np.sum(s*s,axis=(1,2))/(4*G)))
    ratio=rho/2700;volume=float(np.sum(m*BULK/2700*(np.log(ratio)+1/ratio-1)))
    close(r['kinetic_J'],kinetic,1e-25,'actual kinetic energy');close(r['Es_diagnostic_J'],shear+volume,1e-18,'old Es diagnostic')
    c=r['cumulative'];residual=kinetic-K0+shear+volume+potential-initial_u+c['plastic_work_J']+c['friction_raw_J']-c['constraint_work_J']-c['wall_work_J']
    close(r['energy_diagnostic_residual_J'],residual,1e-18,'global diagnostic subtracts initial K')
    warning=abs(residual)>1e-15+.01*(K0+initial_u+abs(c['constraint_work_J'])+abs(c['wall_work_J'])+abs(c['friction_raw_J']))
    require(r['energy_status']==('ENERGY_WARNING' if warning else 'DIAGNOSTIC_WITHIN_SCREEN_NOT_CERTIFIED_FREE_ENERGY'),'raw diagnostic energy warning')
    return dict(x=x,v=v,m=m,rho=rho,s=s,ep=ep,K=kinetic,U=potential,U0=initial_u,P=np.sum(m[:,None]*v,axis=0),
        L=np.sum(m[:,None]*np.cross(x,v),axis=0),normal=normal,control=control,residual=residual,warning=warning)


def audit_trace(trace,index,pid,x,vb,va,force,free,prescribed,dt,step):
    require(trace['particle_index']==index and trace['particle_id']==pid and trace['step']==step,'trace identity')
    require(trace['stage']=='MIDPOINT_FROZEN_GEOMETRY_KICK','trace stage')
    close(trace['time_mid_s'],(step-.5)*dt,1e-20,'trace midpoint')
    close(trace['position_m'],x,2e-17,'trace midpoint geometry')
    close(trace['material_velocity_before_m_s'],vb,1e-12,'trace before velocity')
    close(trace['internal_force_N'],force,1e-12,'trace true internal force')
    close(trace['free_axes'],free,0.,'z-only actuator mobility');close(trace['prescribed_velocity'],prescribed,0.,'actuator command')
    close(trace['mass_kg'],MASS,1e-21,'trace mass')
    n=normal_reference(np.asarray([x]));fn=float(n['force'][0]);faces=trace['contact_rows']
    require(len(faces)==int(n['active'][0]),'actual base face inventory')
    for face in faces:
        require(face['particle_id']==pid and face['feature_id']==0 and face['contact_id']==str(pid)+':0','base contact feature identity')
        close(face['normal'],[0.,0.,1.],0.,'positive-z base normal')
        close(face['application_point'],x,2e-17,'centre force location')
        close(face['contact_point'],[x[0],x[1],-.01],2e-17,'actual base plane')
        close(face['gap'],n['gap'][0],2e-17,'base gap');close(face['penetration'],n['overlap'][0],2e-17,'base overlap')
        close(face['normal_force_N'],fn,1e-12,'gap-derived base force')
        for key in ('wall_velocity','wall_velocity_at_face','translation_velocity','omega'):close(face[key],np.zeros(3),0.,'fixed base '+key)
    result=trace['result'];require(bool(result['intervals']),'accepted kick intervals')
    close(result['constraint_projector'],np.diag(free),0.,'constraint projection')
    close(result['velocity'],va,1e-12,'kick endpoint')
    values=dict(internal=0.,normal=0.,actuator=0.,kinetic=0.);base=np.zeros(3);actuator=np.zeros(3);workv=np.zeros(3)
    elapsed=0.;previous=vb
    for sub in result['intervals']:
        a=audit_interval(sub,previous,force,fn,free,prescribed,MASS,dt,elapsed)
        previous=a['v_after'];elapsed+=a['h'];workv+=a['work_velocity'];base+=a['normal_impulse'];actuator+=a['actuator_impulse']
        for key in values:values[key]+=a[key]
    close(elapsed,dt,1e-20,'whole kick duration');close(previous,va,1e-12,'accepted kick continuity')
    close(result['constraint_impulse'],actuator,1e-17,'summed actuator impulse')
    close(trace['actuator_force_on_material_N'],actuator/dt,1e-10,'reaction force is true impulse/dt')
    close(np.asarray(result['normal_impulse']).reshape(len(faces),3).sum(axis=0),base,1e-18,'summed normal impulse')
    close(np.sum(result['friction_heat_J']),0.,0.,'no fake friction source')
    for target,source in [('other_work_J','internal'),('constraint_work_J','actuator'),('delta_kinetic_J','kinetic')]:
        close(result[target],values[source],1e-18,'summed '+target)
    close(result['wall_work_J'],0.,0.,'stationary base does no wall work')
    return dict(**values,base_impulse=base,actuator_impulse=actuator,workv=workv/dt)


SOURCES=OLD_SOURCES|{'scripts/probe_axial_feed.py','src/afsd_backend/pysph/normal_penalty.py',
    'src/afsd_backend/pysph/axial_actuator.py','configs/benchmarks/axial_feed.json','docs/AXIAL_FEED_CONTRACT.md'}


def audit_vtk(path,r,original):
    vtk=read_point_cloud(path);scalar=vtk['scalars'];step=r['step']
    require(CLASSIFICATION in vtk['title'] and 'NO_THERMAL_SOLUTION' in vtk['title'],'VTK classification')
    names={'particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','boundary_role','as_rs_tag',
        'mass_kg','density_kg_per_m3','pressure_Pa','equivalent_stress_Pa','equivalent_plastic_strain','smoothing_length_m',
        'signed_gap_m','normal_signed_penetration_m','normal_penetration_m','normal_force_N','normal_potential_J','contact_active',
        'actuator_eligible','actuator_active',*('stress_'+a+'_Pa' for a in ('xx','yy','zz','xy','xz','yz'))}
    if step:names|={'equivalent_deviatoric_strain_rate_per_s','equivalent_plastic_strain_rate_per_s'}
    require(set(scalar)==names and set(vtk['vectors'])=={'velocity_m_per_s','initial_position_m'},'VTK solved fields inventory')
    close(vtk['points'],r['position_m'],0.,'VTK position');close(vtk['vectors']['velocity_m_per_s'],r['velocity_m_s'],0.,'VTK velocity')
    close(vtk['vectors']['initial_position_m'],original,0.,'VTK persistent initial position')
    for key in ('particle_id','tracer_type','boundary_role','as_rs_tag'):require(scalar[key]==r[key],'VTK '+key)
    fields=dict(mass_kg='mass_kg',density_kg_per_m3='density_kg_m3',pressure_Pa='pressure_Pa',equivalent_stress_Pa='equivalent_stress_Pa',equivalent_plastic_strain='ep')
    if step:fields.update(equivalent_deviatoric_strain_rate_per_s='strain_rate_per_s',equivalent_plastic_strain_rate_per_s='ep_rate_per_s')
    for target,source in fields.items():close(scalar[target],r[source],0.,'VTK '+target)
    for a,b,key in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        close(scalar['stress_'+key+'_Pa'],np.asarray(r['stress_Pa'])[:,a,b],0.,'VTK stress '+key)
    n=normal_reference(r['position_m']);a=actuator_reference(original[:,2],r['time_s'])
    for key,value in [('signed_gap_m',n['gap']),('normal_signed_penetration_m',n['signed']),('normal_penetration_m',n['overlap']),
        ('normal_force_N',n['force']),('normal_potential_J',n['potential']),('contact_active',n['active'].astype(int)),
        ('actuator_eligible',a['eligible'].astype(int)),('actuator_active',a['active'].astype(int))]:
        close(scalar[key],value,1e-12 if key=='normal_force_N' else 1e-17,'VTK '+key)
    close(scalar['material_id'],np.ones(64),0.,'same material for all tags');close(scalar['component_id'],np.ones(64),0.,'deformable local bar component')
    close(scalar['actuator'],a['eligible'].astype(int),0.,'initial eligibility metadata')
    close(scalar['wholebar_tag'],r['tracer_type'],0.,'whole local probe tags not paper band')
    close(scalar['smoothing_length_m'],np.full(64,.00065),1e-18,'fixed smoothing length')


def audit_case(directory):
    directory=Path(directory);m=load(directory/'manifest.json')
    require(m['status']=='COMPLETE_AXIAL_FEED_COMPONENT_NOT_AFSD','full axial trajectory not complete')
    require(m['production_authorized'] is False and m['mode']=='FULL_REGISTERED_AXIAL_COMPONENT','component scope')
    cfg=m['configuration'];dt=configuration(cfg);count=round(END/dt)
    require(canonical_hash(cfg)==m['configuration_sha256'],'configuration hash')
    require(set(m['source_sha256'])==SOURCES,'complete source inventory')
    for name,digest in m['source_sha256'].items():require(file_hash(ROOT/name)==digest,'source drift '+name)
    files={'history.json','steps.jsonl','last_accepted.npz','release_accepted.npz',*(f'frame_{k:04d}.vtk' for k in range(count+1))}
    actual={p.name for p in directory.iterdir() if p.is_file() and p.name!='manifest.json'}
    require(actual==set(m['artifact_sha256'])==files,'complete artifact inventory')
    for name,digest in m['artifact_sha256'].items():require(file_hash(directory/name)==digest,'artifact drift '+name)
    require(m['accepted_steps']==m['target_steps']==count,'registered steps')
    close(m['accepted_time_s'],END,1e-20,'common full endpoint')
    require(0<m['wall_seconds']<=m['maximum_wall_seconds']<=170,'wall budget')
    h=load(directory/'history.json');journal=[json.loads(line) for line in (directory/'steps.jsonl').read_text(encoding='utf-8').splitlines()]
    require(len(h)==count+1 and len(journal)==count,'history/journal completeness')
    original=literal_geometry();ids=list(range(1,65));eligible=original[:,2]>-.0085;roles=(2*eligible).tolist()
    v0=np.zeros((64,3));v0[eligible,2]=-SPEED;initial_L=np.sum(MASS*np.cross(original,v0),axis=0)
    state=history_record(h[0],original,ids,0,dt)
    for actual,expected,name in [(state['x'],original,'geometry'),(state['v'],v0,'prepared velocity'),
        (state['s'],np.zeros((64,3,3)),'stress'),(state['rho'],np.full(64,2700.),'density'),(state['ep'],np.zeros(64),'plastic strain')]:
        close(actual,expected,0.,'initial '+name)
    require(h[0]['rate_field_status']=='NOT_EVALUATED' and h[0]['strain_rate_per_s'] is None,'unevaluated initial derivative not fabricated')
    sums={key:0. for key in h[0]['cumulative']};require(all(value==0 for value in h[0]['cumulative'].values()),'initial ledger not initial preparation work')
    audit_vtk(directory/'frame_0000.vtk',h[0],original)
    maxima=dict(condition=0.,power_W=0.,internal_force_N=0.,free_speed_m_s=0.,energy_J=0.,xz_Pa=0.,yz_Pa=0.)
    warnings=0;first_contacts={};events=[];active_any=False
    for step,(row,r) in enumerate(zip(journal,h[1:]),1):
        require(row['step']==step and row['particle_id']==ids and row['boundary_role']==roles,'journal identity')
        require(row['stage']=='ACCEPTED_DRIFT_KICK_DRIFT_COMPLETE','accepted journal only')
        close(row['dt_s'],dt,0.,'journal dt');close(row['time_mid_s'],(step-.5)*dt,1e-20,'journal midpoint time')
        current=history_record(r,original,ids,step,dt)
        xm=array(row['x_mid_m'],(64,3),'midpoint geometry');vb=array(row['v_before_m_s'],(64,3),'kick before');va=array(row['v_after_m_s'],(64,3),'kick after')
        close(vb,state['v'],1e-12,'unmodified pre-kick velocity');close(va,current['v'],1e-12,'accepted post-kick velocity')
        close(xm,state['x']+.5*dt*vb,2e-17,'first drift');close(current['x'],xm+.5*dt*va,2e-17,'second drift')
        close(row['x_before_m'],state['x'],0.,'saved pre-step position');close(row['x_after_m'],current['x'],0.,'saved accepted position')
        for name,x,t in [('before',state['x'],(step-1)*dt),('midpoint',xm,(step-.5)*dt),('after',current['x'],step*dt)]:
            audit_normal_snapshot(row['normal_'+name],x,ids,t);audit_control(row['control_'+name],original[:,2],t)
        command=actuator_reference(original[:,2],(step-.5)*dt)
        rho=array(row['rho_before_kg_m3'],(64,),'operator density');close(rho,state['rho'],1e-12,'density chronology')
        close(row['force_density_kg_m3'],rho,0.,'force/operator paired volumes');close(row['mass_kg'],current['m'],0.,'journal mass')
        i,j=np.asarray(row['i'],int),np.asarray(row['j'],int);g=np.asarray(row['g']);B=np.asarray(row['B'])
        sigma=array(row['sigma_mid_Pa'],(64,3,3),'midpoint stress');force=array(row['internal_force_N'],(64,3),'internal force')
        L=reconstruct_gradient(vb,i,j,g,current['m']/rho,B);close(row['L_before_per_s'],L,1e-10,'precontact gradient')
        drho=-rho*np.trace(L,axis1=1,axis2=2);close(row['drho_kg_m3_s'],drho,1e-6,'continuity rate')
        close(current['rho'],rho+dt*drho,1e-10,'density kick');close(row['rho_after_kg_m3'],current['rho'],0.,'saved accepted density')
        close(row['rho_mid_kg_m3'],rho+.5*dt*drho,1e-10,'pressure midpoint density')
        close(row['S_before_Pa'],state['s'],1e-6,'stress chronology');close(row['S_after_Pa'],current['s'],1e-6,'stress accepted state')
        close(sigma,.5*(state['s']+current['s'])-BULK*((rho+.5*dt*drho)/2700-1)[:,None,None]*np.eye(3),1e-6,'force midpoint stress')
        sound=np.sqrt((BULK+4*G/3)/2700);acoustic=.1*.00065/(sound+np.max(np.linalg.norm(vb,axis=1)))
        require(dt<=min(acoustic,.1*np.sqrt(MASS/KN)),'acoustic and penalty timestep bounds')
        traces=row['traces'];require(len(traces)==64,'all particle contact/constraint traces')
        workv=np.zeros((64,3));base=np.zeros((64,3));actuator=np.zeros((64,3));inc=dict.fromkeys(sums,0.)
        for index,trace in enumerate(traces):
            require(trace['actuator_active']==bool(command['active'][index]),'trace active actuator selection')
            a=audit_trace(trace,index,ids[index],xm[index],vb[index],va[index],force[index],
                command['free_axes'][index],command['prescribed_velocity'][index],dt,step)
            workv[index]=a['workv'];base[index]=a['base_impulse'];actuator[index]=a['actuator_impulse']
            for target,source in [('internal_work_J','internal'),('normal_work_J','normal'),('constraint_work_J','actuator'),('delta_kinetic_J','kinetic')]:inc[target]+=a[source]
        close(row['work_velocity_m_s'],workv,1e-12,'actual subinterval work velocity')
        close(row['endpoint_midpoint_velocity_m_s'],.5*(vb+va),1e-12,'separate endpoint midpoint')
        op=audit_operator(xm,workv,current['m'],rho,i,j,g,B,np.asarray(row['L_midpoint_per_s']),force,sigma)
        maxima['condition']=max(maxima['condition'],op['maximum_condition']);maxima['power_W']=max(maxima['power_W'],op['power_residual_W'])
        maxima['internal_force_N']=max(maxima['internal_force_N'],float(np.max(np.linalg.norm(force,axis=1))))
        inc['material_old_work_J']=float(dt*np.einsum('n,nij,nij',current['m']/rho,sigma,L));inc['midpoint_stress_work_J']=dt*op['stress_power_W']
        close(current['ep'],np.zeros(64),0.,'short elastic trajectory only')
        close(row['plastic_work_J'],np.zeros(64),0.,'zero actual plastic work');close(row['return_dissipation_J'],np.zeros(64),0.,'zero actual return dissipation')
        D=.5*(L+L.swapaxes(1,2))-np.trace(L,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
        close(r['strain_rate_per_s'],np.sqrt(2/3*np.sum(D*D,axis=(1,2))),1e-10,'precontact strain-rate field')
        close(r['ep_rate_per_s'],np.zeros(64),0.,'actual plastic strain-rate field')
        inc['normal_relative_work_J']=inc['normal_work_J'];du=current['U']-state['U']
        inc['normal_temporal_residual_J']=du+inc['normal_relative_work_J']
        for key,value in [('potential_change_J',du),('material_work_J',inc['normal_work_J']),('wall_work_J',0.),('temporal_residual_J',inc['normal_temporal_residual_J'])]:
            close(row['normal_energy'][key],value,1e-18,'normal ledger '+key)
        require(row['normal_energy']['status']=='RAW_NOT_HEAT','normal defect not residual heating')
        close(inc['internal_work_J']+inc['midpoint_stress_work_J'],0.,1e-18,'spatial work identity')
        for body,impulse in [('base',base),('actuator',actuator)]:
            total=impulse.sum(axis=0);angular=np.cross(xm,impulse).sum(axis=0)
            close(row['external_impulse_Ns'][body],total,1e-17,'step '+body+' impulse')
            close(row['external_angular_impulse_Nms'][body],angular,1e-19,'step '+body+' angular impulse')
            for k,axis in enumerate('xyz'):
                inc[f'{body}_impulse_{axis}_Ns']=float(total[k]);inc[f'{body}_angular_impulse_{axis}_Nms']=float(angular[k])
        close(current['P']-state['P'],(base+actuator).sum(axis=0),1e-17,'linear momentum increment')
        close(current['L']-state['L'],np.cross(xm,base+actuator).sum(axis=0),1e-19,'angular momentum increment')
        for key,value in inc.items():
            close(row['increments'][key],value,1e-18,'saved increment '+key)
            sums[key]+=value;close(r['cumulative'][key],sums[key],1e-17,'cumulative '+key)
        total_p=np.zeros(3);total_l=np.zeros(3)
        for body in ('base','actuator'):
            dp=np.array([sums[f'{body}_impulse_{axis}_Ns'] for axis in 'xyz']);dl=np.array([sums[f'{body}_angular_impulse_{axis}_Nms'] for axis in 'xyz'])
            close(r[f'cumulative_{body}_impulse_Ns'],dp,1e-17,'reported '+body+' impulse')
            close(r[f'cumulative_{body}_angular_impulse_Nms'],dl,1e-19,'reported '+body+' angular impulse')
            total_p+=dp;total_l+=dl
        close(current['P']-P0,total_p,1e-17,'global momentum subtracts initial P')
        close(current['L']-initial_L,total_l,1e-19,'global angular momentum subtracts initial L')
        close(r['linear_momentum_kg_m_s'],current['P'],1e-22,'actual reported momentum')
        close(r['angular_momentum_kg_m2_s'],current['L'],1e-24,'actual reported angular momentum')
        kinetic=current['K']-K0-sums['constraint_work_J']-sums['internal_work_J']-sums['normal_relative_work_J']
        close(kinetic,0.,1e-18,'global kinetic ledger subtracts initial K');close(r['kinetic_contact_residual_J'],kinetic,1e-18,'saved kinetic ledger')
        require(r['boundary_mass_flux_kg_s']==0. and r['system_mass_status']=='CLOSED_NO_INFINITE_INFLOW','no fictitious mass feed')
        for index,pid in enumerate(ids):
            before=bool(state['normal']['active'][index]);after=bool(current['normal']['active'][index])
            if before!=after:
                event=dict(particle_id=pid,event='CONTACT_ON' if after else 'CONTACT_OFF',bracket_s=[(step-1)*dt,step*dt]);events.append(event)
                if after and str(pid) not in first_contacts:first_contacts[str(pid)]=event['bracket_s']
        active_any=active_any or bool(np.any(current['normal']['active']))
        maxima['free_speed_m_s']=max(maxima['free_speed_m_s'],float(np.max(np.linalg.norm(va[~eligible],axis=1))))
        maxima['xz_Pa']=max(maxima['xz_Pa'],float(np.max(np.abs(current['s'][:,0,2]))));maxima['yz_Pa']=max(maxima['yz_Pa'],float(np.max(np.abs(current['s'][:,1,2]))))
        maxima['energy_J']=max(maxima['energy_J'],abs(current['residual']));warnings+=int(current['warning'])
        audit_vtk(directory/f'frame_{step:04d}.vtk',r,original);state=current
    require(active_any and first_contacts,'actual base contact must occur')
    require(maxima['free_speed_m_s']>1e-8 and maxima['internal_force_N']>1e-6,'force transmission to undriven material')
    require(maxima['xz_Pa']>1 and maxima['yz_Pa']>1,'true 3D off-axis stress response')
    require(m['final_record']==h[-1],'manifest final history identity')
    return dict(status='PASS_AXIAL_FEED_READBACK_NOT_AFSD',classification=CLASSIFICATION,production_authorized=False,
        case=str(directory),dt_s=dt,steps=count,end_time_s=END,frames_checked=count+1,minimum_neighbor_rank=3,
        source_sha256=m['source_sha256'],manifest_sha256=file_hash(directory/'manifest.json'),reader_sha256=file_hash(Path(__file__)),
        independent_operator_reader_sha256=file_hash(ROOT/'scripts/inspect_rotating_deformable_block.py'),
        maxima=maxima,cumulative=sums,initial_kinetic_J=K0,initial_momentum_Ns=P0.tolist(),final_kinetic_J=state['K'],
        final_normal_potential_J=state['U'],final_energy_diagnostic_residual_J=state['residual'],energy_warning_frames=warnings,
        first_contact_brackets_by_id_s=first_contacts,contact_events=events,final_contact_count=int(state['normal']['active'].sum()),
        final_phase=state['control']['phase'],commanded_stroke_m=-SPEED*RELEASE,final_ep_max=float(state['ep'].max()),
        energy_certified=False,continuum_reference='NOT_AVAILABLE; local boundary/adjoint identities not exact continuum trajectory',
        checkpoint='HASH_CHECKED_HERE; actual release restart/continuation tested separately')


def compare_cases(coarse,fine):
    a,b=audit_case(coarse),audit_case(fine)
    require(a['dt_s']==5e-9 and b['dt_s']==2.5e-9 and a['end_time_s']==b['end_time_s']==END,'same registered full endpoint')
    require(a['source_sha256']==b['source_sha256'],'same physics source')
    ca=load(Path(coarse)/'manifest.json')['configuration'];cb=load(Path(fine)/'manifest.json')['configuration']
    require({k:v for k,v in ca.items() if k!='dt_s'}=={k:v for k,v in cb.items() if k!='dt_s'},'same physical problem')
    ha,hb=load(Path(coarse)/'history.json'),load(Path(fine)/'history.json');metrics={}
    for field,limit in dict(position_m=1e-8,velocity_m_s=5e-4,stress_Pa=1e5,density_kg_m3=.1).items():
        differences=[]
        for k,r in enumerate(ha):
            require(r['particle_id']==hb[2*k]['particle_id'],'common frame ID mapping')
            close(r['time_s'],hb[2*k]['time_s'],1e-20,'common physical frame')
            differences.append(float(np.max(np.abs(np.asarray(r[field])-np.asarray(hb[2*k][field])))))
        metrics[field]=dict(maximum_common_difference=max(differences),endpoint_difference=differences[-1],limit=limit,within_screen=max(differences)<=limit)
    for field in ('constraint_work_J','normal_work_J','normal_temporal_residual_J','internal_work_J','friction_raw_J','wall_work_J'):
        differences=[abs(r['cumulative'][field]-hb[2*k]['cumulative'][field]) for k,r in enumerate(ha)]
        metrics[field]=dict(maximum_common_difference=max(differences),endpoint_difference=differences[-1],limit=1e-12,within_screen=max(differences)<=1e-12)
    for field in ('kinetic_J','Es_diagnostic_J','normal_potential_J','energy_diagnostic_residual_J'):
        differences=[abs(r[field]-hb[2*k][field]) for k,r in enumerate(ha)]
        metrics[field]=dict(maximum_common_difference=max(differences),endpoint_difference=differences[-1],limit=1e-12,within_screen=max(differences)<=1e-12)
    contact_comparison={}
    for pid in set(a['first_contact_brackets_by_id_s'])|set(b['first_contact_brackets_by_id_s']):
        aa=a['first_contact_brackets_by_id_s'].get(pid);bb=b['first_contact_brackets_by_id_s'].get(pid)
        contact_comparison[pid]=dict(coarse=aa,fine=bb,upper_sample_difference_s=abs(aa[1]-bb[1]) if aa and bb else None)
    return dict(status='AXIAL_FEED_COMPARISON_COMPLETE_NOT_CONVERGENCE_CERTIFICATE',production_authorized=False,
        screening_status='WITHIN_SCREEN' if all(v['within_screen'] for v in metrics.values()) else 'SENSITIVE_WARNING',
        common_frames=len(ha),common_end_s=END,metrics=metrics,first_contact_comparison_by_id=contact_comparison,cases=[a,b],
        note='Same geometry/boundaries; contact brackets are step-sampled discrete events, not an independent wave arrival reference.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path);parser.add_argument('--coarse',type=Path);parser.add_argument('--fine',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    require(not args.output.exists(),'fresh readback output required')
    if bool(args.case)==bool(args.coarse or args.fine) or (not args.case and not(args.coarse and args.fine)):parser.error('use --case OR --coarse and --fine')
    try:result=audit_case(args.case) if args.case else compare_cases(args.coarse,args.fine);code=0
    except Exception as error:result=dict(status='FAIL_AXIAL_FEED_READBACK',error=str(error),production_authorized=False,reader_sha256=file_hash(Path(__file__)));code=2
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,allow_nan=False));return code


if __name__=='__main__':raise SystemExit(main())
