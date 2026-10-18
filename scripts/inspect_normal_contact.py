"""Independent normal-contact references/readback. Never imports a solver."""
from pathlib import Path
import argparse
import json
import numpy as np
from inspect_rotating_deformable_block import (load, file_hash, close, array, require,
    audit_operator, reconstruct_gradient, canonical_hash)
from afsd_scaffold.vtk_legacy import read_point_cloud

ROOT=Path(__file__).resolve().parents[1]
KN=657500.; MASS=2700*.0005**3; G=26.3e9; BULK=2*G*1.33/(3*.34)
LAYER0=.002-1/KN; U0=16/(2*KN)


def analytic_impact(time, wall_speed):
    """Literal 1kg,4N/m continuous unilateral spring solution."""
    t=np.asarray(time,float); require(np.all(t>=0),'negative reference time')
    active=t<np.pi/2
    delta=np.where(active,.5*np.sin(2*t),0.)
    relative=np.where(active,-np.cos(2*t),1.)
    z=np.where(active,1-delta,1+t-np.pi/2)+wall_speed*t
    impulse=relative+1
    return dict(z=z,v=relative+wall_speed,delta=delta,force=4*delta,potential=2*delta**2,
        normal_impulse=impulse,wall_work=wall_speed*impulse,exit_time=np.pi/2)


def tool_at(t):
    offset=20*max(float(t)-2e-7,0.)
    return LAYER0+offset,20. if t>=2e-7 else 0.


def normal_at(x,t):
    x=np.asarray(x); layer,vz=tool_at(t)
    signed=.00025-(layer-x[:,2]); active=signed>0
    delta=np.where(active,signed,0.)
    return dict(signed=signed,active=active,delta=delta,force=KN*delta,
        potential=.5*KN*delta**2,layer=layer,wall_vz=vz)


def audit_normal_interval(a,vbefore,force,wall,normal_force,free,mass,dt,start,initial):
    h=float(a['dt_s']); require(0<h<=dt,'interval duration')
    close(a['start_s'],start,1e-20,'subinterval chronology')
    vb=array(a['v_before'],(3,),'v_before'); va=array(a['v_after'],(3,),'v_after')
    close(vb,vbefore,1e-12,'subinterval velocity continuity')
    wall=np.asarray(wall,float).reshape(-1,3); nc=len(wall)
    jn=np.asarray(a['normal_impulse'],float).reshape(nc,3)
    jt=np.asarray(a['tangential_impulse'],float).reshape(nc,3)
    reaction=array(a['constraint_impulse'],(3,),'reaction')
    expected=np.tile([0.,0.,-normal_force*h],(nc,1))
    close(jn,expected,1e-18,'gap-derived normal impulse')
    close(jt[:,2],np.zeros(nc),1e-18,'pure tangential impulse')
    close(mass*(va-vb),h*force+jn.sum(axis=0)+jt.sum(axis=0)+reaction,1e-17,'kick momentum')
    close(np.asarray(free)*reaction,np.zeros(3),1e-17,'no hidden free-DOF constraint')
    close((1-np.asarray(free))*va,np.zeros(3),1e-12,'base constraint')
    mid=.5*(vb+va); rel=mid[None,:]-wall; tangent=rel.copy(); tangent[:,2]=0.
    close(np.asarray(a['tangential_relative_midpoint']).reshape(nc,3),tangent,1e-12,'actual tangential midpoint')
    q_each=-np.sum(jt*tangent,axis=1)
    close(a['friction_heat_raw_J'],q_each,1e-18,'raw signed friction')
    saved=np.asarray(a['friction_heat_raw_J']); require(np.all(saved[q_each<0]<0),'friction clipping detected')
    frame=wall[0] if nc else np.asarray(initial)
    close(a['reference_frame_velocity'],frame,1e-12,'wall-frame provenance')
    scale=max(float(np.linalg.norm(initial-frame)),float(np.linalg.norm(wall-frame)) if nc else 0.,
        dt*float(np.linalg.norm(np.asarray(free)*(force+[0,0,-normal_force*nc])/mass)),np.finfo(float).tiny)
    close(a['relative_speed_scale_m_s'],scale,1e-12*max(1.,scale),'independent friction roundoff scale')
    bound=256*np.finfo(float).eps*np.linalg.norm(jt,axis=1)*(scale+np.linalg.norm(vb-frame)+np.linalg.norm(va-frame)+2*np.linalg.norm(wall-frame,axis=1))
    close(a['friction_rounding_bound_J'],bound,max(1e-35,float(np.max(bound,initial=0))*1e-12),'friction rounding bound')
    require(np.all(saved>=-bound),'negative accepted friction')
    if nc:
        cap=.2*normal_force*h; speed=np.linalg.norm(tangent[0])
        require(np.linalg.norm(jt[0])<=cap+1e-18,'friction cap')
        if speed>1e-10: close(jt[0],-cap*tangent[0]/speed,1e-18,'Coulomb sliding direction')
    internal=float(h*force@mid); normal=float(np.sum(jn*mid)); normal_wall=float(np.sum(jn*wall))
    wall_work=float(np.sum((jn+jt)*wall)); constraint=float(reaction@mid)
    kinetic=.5*mass*float(va@va-vb@vb); q=float(q_each.sum())
    residual=kinetic+q-wall_work-constraint-internal-normal+normal_wall
    tol=1e-15+1e-9*(abs(kinetic)+abs(q)+abs(wall_work)+abs(constraint)+abs(internal)+abs(normal)+abs(normal_wall))
    close(residual,0.,tol,'normal-inclusive kinetic work identity')
    for key,expected in [('other_work_J',internal),('constraint_work_J',constraint),('wall_work_J',wall_work),('delta_kinetic_J',kinetic)]:
        close(a[key],expected,tol,'saved '+key)
    close(np.sum(a['normal_work_J']),normal,tol,'saved normal material work')
    return dict(v_after=va,h=h,work_velocity=h*mid,internal=internal,normal=normal,
        normal_wall=normal_wall,wall=wall_work,constraint=constraint,kinetic=kinetic,q=q,
        external_impulse=jn.sum(axis=0)+jt.sum(axis=0)+reaction,contact_impulse=jn.sum(axis=0)+jt.sum(axis=0),
        negative=int(np.sum(q_each<0)),min_q=float(q_each.min()) if nc else None)


def provenance(directory,status):
    directory=Path(directory); m=load(directory/'manifest.json')
    require(m['status']==status,'case not complete')
    require(m['production_authorized'] is False,'production authorization')
    require(canonical_hash(m['configuration'])==m['configuration_sha256'],'config hash')
    for name,expected in m['source_sha256'].items():
        target=(ROOT/name).resolve(); require(ROOT in target.parents,'unsafe source path')
        require(file_hash(target)==expected,'source drift '+name)
    inventory={p.name for p in directory.iterdir() if p.is_file() and p.name!='manifest.json'}
    require(set(m['artifact_sha256'])==inventory,'artifact inventory')
    for name,expected in m['artifact_sha256'].items():
        require(Path(name).name==name,'unsafe artifact path')
        require(file_hash(directory/name)==expected,'artifact drift '+name)
    require(m['wall_seconds']<=m['maximum_wall_seconds'],'wall budget')
    return m


def audit_impact(directory):
    directory=Path(directory); m=provenance(directory,'COMPLETE_MANUFACTURED_NORMAL_IMPACT_NOT_AFSD')
    expected_sources={'scripts/probe_normal_contact.py','configs/benchmarks/normal_contact.json',
        'docs/NORMAL_CONTACT_CONTRACT.md','src/afsd_scaffold/vtk_legacy.py',
        *('src/afsd_backend/pysph/'+name+'.py' for name in ('normal_penalty','rigid_contact_kinematics',
          'geometry','contact','contact_kick','contact_cone','integrator','output')),
        *('third_party/PySPH/pysph/sph/'+name+'.py' for name in ('rigid_body','integrator','integrator_step'))}
    require(set(m['source_sha256'])==expected_sources,'complete impact source inventory')
    require(set(m['artifact_sha256'])=={'history.json','normal_intervals.json','initial.vtk','final.vtk','last_accepted.npz'},'complete impact artifact inventory')
    cfg=m['configuration']; dt=cfg['dt_s']; speed=cfg['wall_translation_velocity_m_s'][2]
    require(dt in (.005,.0025) and speed in (0.,.2),'registered impact case')
    for key,value in [('mass_kg',1.),('radius_m',1.),('stiffness_N_m',4.),('mu',0.),('end_time_s',2.)]:
        require(cfg[key]==value,'literal impact parameter '+key)
    require(cfg['source_configuration']==load(ROOT/'configs/benchmarks/normal_contact.json'),'normal config content')
    literal=dict(mass=1,stiffness=4,radius=1,initial_gap=1,initial_relative_normal_velocity=-1,friction=0,end_time=2)
    require({k:v['value'] for k,v in cfg['source_configuration']['parameters'].items()}==literal,'source impact literals')
    require(cfg['free_axes']==[1,1,1] and cfg['normal']==[0,0,1],'unconstrained positive-z plane')
    close(cfg['initial_position_m'],[0,0,1],0.,'impact initial condition')
    close(cfg['initial_velocity_m_s'],[0,0,speed-1],0.,'impact initial velocity')
    h=load(directory/'history.json'); traces=load(directory/'normal_intervals.json'); count=round(2/dt)
    require(len(h)==count+1 and len(traces)==count and m['accepted_steps']==count,'impact histories complete')
    close(m['accepted_time_s'],2.,1e-14,'impact common endpoint')
    jn_total=0.; work_total=0.; normal_work=0.; defects=0.; maxenergy=0.; maxxerr=0.; maxverr=0.
    exit_bracket=None; ever=False
    for k,r in enumerate(h):
        t=k*dt; close(r['time_s'],t,1e-14,'impact record time')
        require(r['step']==k and r['particle_id']==101,'impact identity')
        require(r['tracer_type']==1 and r['as_rs_tag']==-1,'impact labels')
        x=array(r['position_m'],(3,),'impact position'); v=array(r['velocity_m_s'],(3,),'impact velocity')
        close(x[:2],np.zeros(2),0.,'normal-only position'); close(v[:2],np.zeros(2),0.,'normal-only velocity')
        delta=max(1-(x[2]-speed*t),0.); potential=2*delta*delta; kinetic=.5*float(v@v)
        close(r['gap_m'],x[2]-speed*t,1e-14,'impact gap')
        require(r['active']==(delta>0),'impact unilateral active')
        close(r['normal_force_N'],4*delta,1e-12,'impact force')
        close(r['normal_potential_J'],potential,1e-12,'impact spring potential')
        close(r['kinetic_J'],kinetic,1e-12,'impact kinetic energy')
        ref=analytic_impact(t,speed); maxxerr=max(maxxerr,abs(x[2]-float(ref['z']))); maxverr=max(maxverr,abs(v[2]-float(ref['v'])))
        if k:
            row=traces[k-1]; prior=h[k-1]; vb=np.asarray(prior['velocity_m_s']); x0=np.asarray(prior['position_m'])
            xm=array(row['position_mid_m'],(3,),'impact midpoint'); tm=(k-.5)*dt
            close(xm,x0+.5*dt*vb,1e-13,'impact first drift'); close(x,xm+.5*dt*v,1e-13,'impact second drift')
            force=4*max(1-(xm[2]-speed*tm),0.)
            close(row['normal_state']['normal_force_N'],force,1e-12,'midpoint spring force')
            close(row['wall_velocity_m_s'],[0,0,speed],0.,'moving wall speed')
            elapsed=0.; vv=vb.copy(); this_normal=0.; this_wall=0.; impulse=0.
            for sub in row['result']['intervals']:
                ds=sub['dt_s']; close(sub['start_s'],elapsed,1e-14,'impact interval time')
                require(ds>0 and ds<=dt,'impact subinterval duration')
                before=np.asarray(sub['v_before']); after=np.asarray(sub['v_after']); mid=.5*(before+after)
                close(before,vv,1e-12,'impact kick continuity')
                jn=np.asarray(sub['normal_impulse']).reshape(1,3); jt=np.asarray(sub['tangential_impulse']).reshape(1,3)
                close(jn,[[0,0,force*ds]],1e-12,'spring normal impulse')
                close(jt,np.zeros((1,3)),0.,'zero-mu tangent impulse')
                close(sub['friction_heat_raw_J'],[0.],0.,'zero-mu raw friction')
                close(sub['constraint_impulse'],np.zeros(3),1e-12,'no normal fixture constraint')
                close(after-before,jn[0],1e-12,'mass-one normal momentum')
                wn=float(jn[0]@mid); ww=float(jn[0,2]*speed)
                close(np.sum(sub['normal_work_J']),wn,1e-12,'normal material work')
                close(sub['wall_work_J'],ww,1e-12,'normal moving-wall work')
                close(.5*(after@after-before@before),wn,1e-12,'normal kinetic identity')
                impulse+=jn[0,2]; this_normal+=wn; this_wall+=ww; vv=after; elapsed+=ds
            close(elapsed,dt,1e-14,'whole normal kick duration'); close(vv,v,1e-12,'impact endpoint velocity')
            defect=potential-prior['normal_potential_J']+this_normal-this_wall
            close(row['normal_potential_temporal_defect_J'],defect,1e-12,'raw normal temporal defect')
            jn_total+=impulse; work_total+=this_wall; normal_work+=this_normal; defects+=defect
        c=r['cumulative']
        for key,value in [('normal_impulse_N_s',jn_total),('wall_work_J',work_total),('normal_work_J',normal_work),('normal_potential_temporal_defect_J',defects)]:
            close(c[key],value,1e-11,'impact cumulative '+key)
        close(c['friction_raw_J'],0.,0.,'no manufactured friction heat')
        residual=kinetic+potential-.5*(speed-1)**2-work_total
        close(r['global_energy_residual_J'],residual,1e-12,'physical spring energy defect')
        maxenergy=max(maxenergy,abs(residual))
        if r['active']: ever=True
        elif ever and exit_bracket is None: exit_bracket=[(k-1)*dt,t]
    require(exit_bracket is not None and abs(exit_bracket[1]-np.pi/2)<=2*dt,'sampled spring separation')
    ref=analytic_impact(2.,speed); endpointx=abs(h[-1]['position_m'][2]-float(ref['z'])); endpointv=abs(h[-1]['velocity_m_s'][2]-float(ref['v']))
    require(endpointx<=5e-5 and endpointv<=2e-5 and maxenergy<=1e-4,'analytic impact tolerance')
    close(jn_total,2.,2e-5,'analytic final normal impulse'); close(work_total,2*speed,1e-5,'analytic moving wall work')
    for name,r in [('initial.vtk',h[0]),('final.vtk',h[-1])]:
        vtk=read_point_cloud(directory/name)
        close(vtk['points'],[r['position_m']],0.,'impact VTK position')
        close(vtk['vectors']['velocity_m_per_s'],[r['velocity_m_s']],0.,'impact VTK velocity')
        require(set(vtk['scalars'])=={'particle_id','tracer_type','as_rs_tag','mass_kg'},'impact unsolved fields not fabricated')
        for field,value in [('particle_id',101),('tracer_type',1),('as_rs_tag',-1),('mass_kg',1.)]:
            close(vtk['scalars'][field],[value],0.,'impact VTK '+field)
    require(m['final_record']==h[-1],'impact final state identity')
    return dict(status='PASS_NORMAL_IMPACT_READBACK_NOT_AFSD',production_authorized=False,
        dt_s=dt,wall_speed_m_s=speed,end_time_s=2.,steps=count,endpoint_position_error_m=endpointx,
        endpoint_velocity_error_m_s=endpointv,max_history_position_error_m=maxxerr,max_history_velocity_error_m_s=maxverr,
        max_physical_energy_defect_J=maxenergy,final_physical_energy_defect_J=h[-1]['global_energy_residual_J'],
        separation_time_bracket_s=exit_bracket,analytic_separation_s=float(np.pi/2),normal_impulse_Ns=jn_total,
        wall_work_J=work_total,normal_work_J=normal_work,raw_normal_temporal_defect_J=defects,
        wall_seconds=m['wall_seconds'],particle_array_bytes=m['particle_array_bytes'],
        manifest_sha256=file_hash(directory/'manifest.json'),source_sha256=m['source_sha256'],reader_sha256=file_hash(Path(__file__)))


def compare_impacts(coarse,fine):
    a,b=audit_impact(coarse),audit_impact(fine)
    require(a['dt_s']==.005 and b['dt_s']==.0025,'registered impact dt pair')
    require(a['wall_speed_m_s']==b['wall_speed_m_s'],'same translating wall')
    require(a['source_sha256']==b['source_sha256'],'same impact source')
    ca=load(Path(coarse)/'manifest.json')['configuration']; cb=load(Path(fine)/'manifest.json')['configuration']
    require({k:v for k,v in ca.items() if k!='dt_s'}=={k:v for k,v in cb.items() if k!='dt_s'},'same impact physics')
    ha,hb=load(Path(coarse)/'history.json'),load(Path(fine)/'history.json')
    differences={field:max(float(np.max(np.abs(np.asarray(row[field])-np.asarray(hb[2*k][field]))))
        for k,row in enumerate(ha)) for field in ('position_m','velocity_m_s','normal_potential_J')}
    return dict(status='PASS_ANALYTIC_IMPACT_PAIR_NOT_AFSD',production_authorized=False,
        common_end_s=2.,common_frames=len(ha),cases=[a,b],common_history_max_absolute=differences,
        note='Each run independently meets preregistered analytic limits; two dt values do not certify an asymptotic order.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--impact',type=Path); parser.add_argument('--coarse',type=Path); parser.add_argument('--fine',type=Path)
    parser.add_argument('--output',type=Path,required=True); args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('fresh readback output required')
    if bool(args.impact)==bool(args.coarse or args.fine) or (not args.impact and not(args.coarse and args.fine)):
        parser.error('use --impact OR --coarse and --fine')
    try:
        result=audit_impact(args.impact) if args.impact else compare_impacts(args.coarse,args.fine); code=0
    except Exception as error:
        result=dict(status='FAIL_NORMAL_READBACK',error=str(error),production_authorized=False,reader_sha256=file_hash(Path(__file__))); code=2
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,allow_nan=False)); return code


if __name__=='__main__': raise SystemExit(main())
