"""Literal slow-load manufactured inputs ONLY, shared with independent oracle."""
import copy
from reference_hht_block import fixture

CLASSIFICATION='LOCAL_SLOW_LOAD_FULL_STARTUP_FEASIBILITY_NOT_AFSD'
SCOUT_END_S=20e-6
DESIRED_END_S=1e-3
REFERENCE_STEPS_S=(7.8125e-9,3.90625e-9)
SNAPSHOT_STEP_S=62.5e-9


def inputs():
    data=copy.deepcopy(fixture())
    data['classification']=CLASSIFICATION
    data['parameters'].update(end_time=SCOUT_END_S,desired_end_time_s=DESIRED_END_S,
        normal_frequency_Hz=5.,wall_tangential_program='CONSTANT_PLUS_Z',
        maximum_translation_m=250e-6,maximum_pair_relative_change=.005)
    data['program_source']='MANUFACTURED_SLOW_LOADING_NOT_AFSD_TOOL; SEE PROCESS_SLOW_CONTRACT'
    return data


def schedule(refinement=0):
    if type(refinement) is not int or refinement not in (0,1,2):
        raise ValueError('registered schedules S0/S1/S2 only')
    result=[];begin=0.
    for end,base_count in ((2e-6,8),(4e-6,4),(20e-6,16)):
        count=base_count*2**refinement;step=(end-begin)/count
        for index in range(count):result.append((begin+index*step,begin+(index+1)*step))
        begin=end
    return result
