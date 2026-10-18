"""Independent manufactured control-volume identities, not coupled solves."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import inlet_protocol as inventory
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph import interval_acceptance as crossing
import inlet_acceptance_v2 as acceptance
from inlet_source_work import evaluate_source_work,source_work_record,PORT_CONTRACT


def fixture(reference=300.,birth_velocity=(0.,0.,0.),raw_defect=0.,unload=False):
    p=dict(dx=.1,h_ratio=1.3,rho0=2.,width=1.,marker_low=0.,marker_high=1.,bar_length=1.,
        actuator_depth=.1,temperature=310.,traverse=0.,feed=.1,shear=3.,poisson=.25,cp=10.,
        plastic_heat_fraction=.8,heat_partition=.6)
    parent=make_array([[0.,0.,.1],[.1,.1,.2]],np.ones(2,int),p);parent.m[:]=1.
    if unload:parent.S[:9]=np.diag([2.,-2.,0.]).ravel()
    inc=make_array([[0.,0.,.5]],np.ones(1,int),p);inc.particle_id[:]=0;inc.m[:]=1.;inc.temperature[:]=320.
    inc.u[:],inc.v[:],inc.w[:]=birth_velocity
    fields={key:dict(unit=inventory.REQUIRED.get(key,'1'),source='MANUFACTURED_INVENTORY_INITIAL_CONDITION',
        status='INITIAL_CONDITION' if key in inventory.REQUIRED else 'UNRESOLVED') for key in inc.properties}
    spec=dict(event_id='manufactured-birth-1',time_s=0.,source='MANUFACTURED_STRESS_FREE_RESERVOIR',fields=fields,
        heat_capacity=dict(value=10.,unit='J/(kg K)',source='MANUFACTURED_CONSTANT',status='DECLARED_ASSUMPTION'),
        reference_temperature=dict(value=reference,unit='K',source='MANUFACTURED_REFERENCE',status='DECLARED_ASSUMPTION'))
    registry=inventory.initialize_registry(parent,'MANUFACTURED_COMPONENT_ONLY',acknowledgement=inventory.ACKNOWLEDGEMENT)
    prepared,_,event=inventory.propose_insertion(parent,registry,inc,spec)
    prepared.actuator[-1]=1
    actual=prepared.extract_particles(np.arange(3,dtype=np.uint32))
    displacement=np.zeros((3,3));displacement[-1,2]=-.01
    actual.x[:],actual.y[:],actual.z[:]=(np.column_stack((prepared.x,prepared.y,prepared.z))+displacement).T
    actual.u[:],actual.v[:],actual.w[:]=(displacement/.1).T
    free=np.ones((3,3),bool);free[0]=False;free[-1,2]=False
    if unload:
        actual.S[:]=0.;dU=-1/3;plastic=friction=normal=tangent=tool=contact_in=0.
    else:
        actual.S[:9]=np.diag([2.,-2.,0.]).ravel();dU=1/3
        plastic=.3;friction=.1;normal=.01;tangent=.02;tool=.04;contact_in=-.05
    external_heat=.8*plastic+.6*friction+contact_in
    actual.temperature[0]+=external_heat/10.
    measured=float(np.sum(10.*(actual.temperature-prepared.temperature)))
    Wact=dU+normal+tangent+plastic+friction-tool+raw_defect
    force=np.zeros((3,3));force[-1,2]=Wact/.01
    contact={key:np.zeros((3,8)) for key in ('normal_potential0_J','normal_potential1_J',
        'tangential_potential0_J','tangential_potential1_J','wall_work_J','friction_dissipation_J')}
    contact['normal_potential1_J'][0,0]=normal;contact['tangential_potential1_J'][0,0]=tangent
    contact['wall_work_J'][0,0]=tool;contact['friction_dissipation_J'][0,0]=friction
    wp=np.array([plastic/.5,0.,0.])
    thermal=dict(plastic_work_J=plastic,plastic_heat_J=.8*plastic,plastic_unthermalized_J=.2*plastic,
        friction_raw_J=friction,friction_material_heat_J=.6*friction,friction_wall_heat_J=.4*friction,
        measured_thermal_increment_J=measured,wall_heat_out_J=.4*friction-contact_in,
        contact_in_J=contact_in,internal_conduction_J=0.,mc_delta_T_residual_J=measured-external_heat)
    result=dict(displacement=displacement,position=np.column_stack((actual.x,actual.y,actual.z)),
        force=force,density=actual.rho.copy(),stress=actual.S.reshape(3,3,3).copy(),temperature=actual.temperature.copy(),
        volume=np.full(3,.5),dt_s=.1,material=(actual.S.reshape(3,3,3).copy(),actual.ep.copy(),wp,np.zeros(3),np.zeros(3)),
        contact=contact,thermal=dict(ledger=thermal))
    energy=dict(material_storage_increment_J=dU,plastic_work_J=plastic,normal_storage_increment_J=normal,
        tangential_storage_increment_J=tangent,friction_sliding_J=friction,prescribed_boundary_work_J=Wact,
        tool_work_J=tool,global_raw_residual_J=raw_defect,prescribed_reaction_N=-force*(~free))
    port=dict(contract=PORT_CONTRACT,surface='EXTERIOR_BIRTH_BOUNDARY_NOT_INTERNAL_CONTROL_PLANE',
        internal_interface_work='NOT_ADDED_INTERNAL_PAIR_WORK',traction_Pa=[0.,0.,0.],
        source='MANUFACTURED_ZERO_TRACTION_EXTERIOR_BOUNDARY')
    args=dict(parameters=p,free=free,result=result,mechanical_energy=energy,parent_contact_storage_J=np.zeros((2,2)),external_port=port)
    return parent,prepared,actual,event,args


def evaluate(data):return evaluate_source_work(*data[:4],**data[4])


class InletSourceWorkTests(unittest.TestCase):
    def test_independent_heated_elastic_contact_control_volume_identity(self):
        report=evaluate(fixture());energy=report['open_energy'];terms=report['computed_terms']
        # m=1 kg, rho=2, G=3 Pa, S=diag(2,-2,0): U=(1/2)*8/12=1/3 J.
        # Birth carries 1*10*(320-300)=200 J, no elastic energy.
        self.assertEqual(report['source_increment']['mass_kg'],1.)
        self.assertEqual(report['source_increment']['sensible_heat_J'],200.)
        self.assertEqual(terms['elastic_convected_energy']['value'],0.)
        self.assertEqual(terms['pressure_flow_work']['value'],0.)
        self.assertAlmostEqual(terms['mechanical_reaction_work']['value'],1/3+.01+.02+.3+.1-.04)
        self.assertAlmostEqual(energy['material_storage_after_J'],1/3)
        self.assertAlmostEqual(energy['total_storage_increment_J'],200.+.25+1/3+.03)
        self.assertAlmostEqual(energy['wall_heat_out_J'],.09)
        self.assertAlmostEqual(energy['plastic_unthermalized_J'],.06)
        self.assertAlmostEqual(energy['open_raw_residual_J'],0.,delta=2e-12)
        self.assertEqual(report['status'],'COMPUTED_SOURCE_WORK_NOT_PROCESS_ACCEPTANCE')

    def test_negative_reaction_work_unloading_is_not_clipped(self):
        report=evaluate(fixture(unload=True))
        self.assertAlmostEqual(report['computed_terms']['mechanical_reaction_work']['value'],-1/3)
        self.assertAlmostEqual(report['open_energy']['open_raw_residual_J'],0.,delta=2e-12)

    def test_raw_defect_survives_without_residual_heat(self):
        base=evaluate(fixture());bad=evaluate(fixture(raw_defect=.012))
        self.assertAlmostEqual(bad['open_energy']['mechanical_raw_residual_J'],.012)
        self.assertAlmostEqual(bad['open_energy']['open_raw_residual_J'],.012,delta=2e-12)
        self.assertEqual(base['open_energy']['sensible_storage_after_J'],bad['open_energy']['sensible_storage_after_J'])

    def test_reference_temperature_mass_source_shift_cancels(self):
        a=evaluate(fixture(reference=300.));b=evaluate(fixture(reference=290.))
        self.assertEqual(b['source_increment']['sensible_heat_J']-a['source_increment']['sensible_heat_J'],100.)
        self.assertAlmostEqual(b['open_energy']['total_storage_increment_J']-a['open_energy']['total_storage_increment_J'],100.)
        self.assertAlmostEqual(a['open_energy']['open_raw_residual_J'],b['open_energy']['open_raw_residual_J'],delta=2e-12)
        self.assertEqual(b['thermal_inventory_reference']['reference_temperature_K'],290.)

    def test_incoming_momentum_and_kinetic_energy_do_not_enter_qs_heat(self):
        quiet=evaluate(fixture());moving=evaluate(fixture(birth_velocity=(1.,2.,3.)))
        self.assertEqual(moving['source_increment']['momentum_kg_m_s'],[1.,2.,3.])
        self.assertEqual(moving['source_increment']['kinetic_energy_J'],7.)
        self.assertFalse(moving['quasistatic_kinetic_diagnostic']['included_in_quasistatic_energy'])
        self.assertFalse(moving['quasistatic_kinetic_diagnostic']['converted_to_heat'])
        self.assertEqual(quiet['open_energy'],moving['open_energy'])

    def test_nonzero_traction_stress_density_and_birth_contact_fail_closed(self):
        for mode in ('traction','stress','rho','history','contact','internal-work'):
            data=fixture();prepared=data[1];args=data[4]
            if mode=='traction':args['external_port']['traction_Pa'][2]=1.
            elif mode=='stress':prepared.S[-9]=.01
            elif mode=='rho':prepared.rho[-1]=2.01
            elif mode=='history':prepared.ep[-1]=.1
            elif mode=='contact':args['result']['contact']['normal_potential0_J'][-1,0]=1e-8
            else:args['external_port']['internal_interface_work']='ADD_INTERNAL_CUT_PRESSURE_WORK'
            with self.subTest(mode=mode),self.assertRaises(ValueError):evaluate(data)

    def test_old_inventory_mutation_extra_work_and_false_ledgers_reject(self):
        for mode in ('old-field','extra-boundary','thermal','mechanical','source','actual'):
            data=fixture();args=data[4]
            if mode=='old-field':data[1].temperature[0]+=1.
            elif mode=='extra-boundary':args['result']['displacement'][0,0]=.01
            elif mode=='thermal':args['result']['thermal']['ledger']['wall_heat_out_J']+=.01
            elif mode=='mechanical':args['mechanical_energy']['prescribed_boundary_work_J']*=2.
            elif mode=='source':data[3]['source_increment']['mass_kg']+=1.
            else:data[2].temperature[0]+=1.
            with self.subTest(mode=mode),self.assertRaises(ValueError):evaluate(data)

    def test_actual_sorting_and_input_purity(self):
        data=fixture();a=evaluate(data)
        before=[{key:prop.get_npy_array().copy() for key,prop in pa.properties.items()} for pa in data[:3]]
        changed=list(data);changed[2]=data[2].extract_particles(np.array([2,0,1],np.uint32))
        b=evaluate(tuple(changed));self.assertEqual(a['open_energy'],b['open_energy'])
        for pa,saved in zip(data[:3],before):
            for key,value in saved.items():np.testing.assert_array_equal(pa.properties[key].get_npy_array(),value)

    def test_artifact_bound_v2_source_work_is_computed_and_tamper_rejected(self):
        report=evaluate(fixture())
        with tempfile.TemporaryDirectory(prefix='inlet_source_work_component_') as directory:
            path=Path(directory)/'manufactured_source_work.json';path.write_text(json.dumps(report),encoding='utf-8')
            record=source_work_record(report,path)
            self.assertEqual(record['mechanical_reaction_work']['value'],report['computed_terms']['mechanical_reaction_work']['value'])
            self.assertEqual(len(record['mechanical_reaction_work']['evidence_sha256']),64)
            corrupted=copy.deepcopy(report);corrupted['computed_terms']['pressure_flow_work']['value']=1.
            path.write_text(json.dumps(corrupted),encoding='utf-8')
            with self.assertRaises(ValueError):source_work_record(report,path)
            with self.assertRaises(ValueError):source_work_record(corrupted,path)

    def test_computed_record_binds_to_diagnostic_v2_transaction(self):
        # Entire coupled map/receipt is manufactured here. This checks that a
        # genuinely computed source record can replace the old fixture values.
        parent,prepared,actual,event,args=fixture()
        tool_pose=lambda _:dict(origin=np.zeros(3),rotation=np.eye(3))
        observation=crossing.initial_bore_crossing_state(np.column_stack((parent.x,parent.y,parent.z)),
            parent.particle_id,parent.component_id,parent.m,0.,tool_pose(0.),1.)
        registry=inventory.initialize_registry(parent,'MANUFACTURED_ADAPTER_REFERENCE',acknowledgement=inventory.ACKNOWLEDGEMENT)
        state=acceptance.initialize_acceptance_state(parent,registry,dict(time_s=0.,contact_history=np.zeros((2,8,3)),
            free=args['free'][:2],bore_observation=observation),
            dict(frame='WORLD_Z',plane_z_m=.4,unit='m',source='MANUFACTURED_CONTROL_PLANE'),'MANUFACTURED_ADAPTER')
        inc=prepared.extract_particles(np.array([2],np.uint32));inc.particle_id[:]=0
        spec=dict(event_id=event['event_id'],time_s=0.,source=event['source'],fields=event['field_declarations'],
            heat_capacity=event['heat_capacity'],reference_temperature=event['reference_temperature'])
        candidate,pending,new_event=acceptance.propose_insertion(parent,state,inc,spec,tool_pose(0.))
        report=evaluate_source_work(parent,candidate,actual,new_event,**args)
        with tempfile.TemporaryDirectory(prefix='inlet_source_v2_adapter_') as directory:
            path=Path(directory)/'manufactured_source_work.json';path.write_text(json.dumps(report),encoding='utf-8')
            terms=source_work_record(report,path);history=np.zeros((3,8,3))
            receipt=dict(scope=acceptance.DIAGNOSTIC_SCOPE,step_id='manufactured-adapter-step',
                binding=acceptance.mechanical_binding(actual,pending,history,args['free'],.1),
                source_work_sha256=acceptance._sha(terms),
                gates={key:dict(status='PASS',evidence_sha256='a'*64) for key in acceptance.GATES})
            final=acceptance.seal_mechanical_acceptance(actual,pending,history,args['free'],.1,receipt,terms,
                lambda _:True,tool_pose)
        self.assertEqual(final['scope_lock'],acceptance.DIAGNOSTIC_SCOPE)
        self.assertFalse(final['process_acceptance_by_this_module'])
        self.assertEqual(final['source_work_events'][0]['work'],terms)

    def test_large_prior_storage_small_increment_consistency_not_false_failure(self):
        data=fixture();parent,prepared,actual,_,args=data
        large=np.sqrt(12.e8);small=np.sqrt(12.e-8)
        old=np.array([np.diag([large,-large,0.]),np.diag([small,-small,0.])])
        parent.S[:]=old.ravel();prepared.S[:18]=old.ravel();actual.S[:18]=old.ravel()
        actual.S[-9:]=np.diag([small,-small,0.]).ravel()
        data[3]['parent_particle_state_sha256']=inventory._array_digest(parent)
        # Reproduce the ORIGINAL core summation order independently, not the
        # evaluator's fsum. Small terms below large-stock ulps expose the gap.
        core_start=np.sum((prepared.m/prepared.rho)*np.sum(prepared.S.reshape(3,3,3)**2,axis=(1,2))/(4.*3.))
        core_end=np.sum((actual.m/actual.rho)*np.sum(actual.S.reshape(3,3,3)**2,axis=(1,2))/(4.*3.))
        increment=float(core_end-core_start)
        energy=args['mechanical_energy'];energy['material_storage_increment_J']=increment
        work=increment+.01+.02+.3+.1-.04
        args['result']['force'][-1,2]=work/.01;energy['prescribed_boundary_work_J']=work
        energy['prescribed_reaction_N']=-args['result']['force']*(~args['free'])
        args['result']['stress']=actual.S.reshape(3,3,3).copy()
        report=evaluate(data)
        rows={row['field']:row for row in report['mechanical_ledger_consistency_roundoff']}
        row=rows['material_storage_increment_J']
        self.assertGreater(abs(row['provided_J']-row['recomputed_J']),1e-10)
        self.assertGreater(row['absolute_storage_terms_scale_J'],1e8)
        self.assertFalse(row['acceptance_threshold'])
        self.assertEqual(energy['global_raw_residual_J'],0.)
        self.assertEqual(report['open_energy']['provided_mechanical_raw_residual_J'],0.)
        energy['material_storage_increment_J']+=.01
        with self.assertRaisesRegex(ValueError,'actual mechanical ledger material_storage_increment_J'):evaluate(data)


if __name__=='__main__':unittest.main(verbosity=2)
