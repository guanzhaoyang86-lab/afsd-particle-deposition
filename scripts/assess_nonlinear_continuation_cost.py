"""Read-only geometric/time-scale feasibility, not a process-time simulation."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np
from pysph.solver.utils import load

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def geometric_gap(z0,layer,guard):
    clearance=np.asarray(z0,float)-float(layer)
    if clearance.ndim!=1 or not clearance.size or not np.isfinite(clearance).all() or guard<0:
        raise ValueError('invalid original bar geometry')
    above=clearance>=0
    if not np.any(above):raise ValueError('no upstream bar cohort; nominal crossing-time estimate undefined')
    return dict(minimum_upstream_original_gap_m=float(clearance[above].min()),
        minimum_absolute_original_gap_m=float(np.abs(clearance).min()),
        initially_on_or_above_plane_count=int(above.sum()),initially_below_plane_count=int((~above).sum()),
        maximum_allowed_displacement_m=float(guard),
        any_center_crossing_permitted_by_current_guard=bool(np.any(np.abs(clearance)<=guard)),
        minimum_absolute_gap_under_guard_m=float(np.abs(clearance).min()-guard),
        proof='All centres start farther than dmax from the plane; norm(x-x_original)<=dmax prohibits crossing in either direction.',
        interpretation='Initially below-plane material is not a newly observed crossing event.')


def assess(case):
    case=Path(case).resolve()
    if not case.is_relative_to(ROOT):raise ValueError('project-local case required')
    m=json.loads((case/'manifest.json').read_text(encoding='utf-8'))
    if m['status']!='COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION' or m['accepted_steps']!=4:
        raise ValueError('registered complete four-stage evidence required')
    p=json.loads((case/'resolved_config.json').read_text(encoding='utf-8'))['parameters']
    source=(ROOT/m['parent_reference']['path']).resolve()
    if not source.is_relative_to(ROOT) or sha(source)!=m['parent_reference']['sha256']:
        raise ValueError('source checkpoint provenance mismatch')
    source_hashes={str((case/'manifest.json').relative_to(ROOT)):sha(case/'manifest.json'),
        str((case/'resolved_config.json').relative_to(ROOT)):sha(case/'resolved_config.json'),
        str(source.relative_to(ROOT)):sha(source)}
    pa=load(str(source))['arrays']['solid'];bar=pa.component_id==1
    original_clearance=pa.z0[bar]-p['layer'];current_clearance=pa.z[bar]-p['layer']
    geometry=geometric_gap(pa.z0[bar],p['layer'],p['maximum_displacement_m'])
    gap=geometry['minimum_upstream_original_gap_m']
    if p['feed']<=0:raise ValueError('positive feed required for nominal actuator-scale estimate')
    nominal=gap/p['feed'];dt=m['settings']['dt_s'];steps=math.ceil(nominal/dt)
    mean_cost=m['wall_seconds']/m['accepted_steps'];projection=steps*mean_cost
    result=dict(status='READ_ONLY_PROCESS_SCALE_GAP_QUANTIFIED_NOT_AUTHORIZATION',production_authorized=False,
        input_sha256=source_hashes,source_sha256=sha(__file__),
        actual=dict(accepted_steps=m['accepted_steps'],physical_increment_s=m['accepted_time_s']-m['source_time_s'],
            end_time_s=m['accepted_time_s'],wall_seconds=m['wall_seconds'],dt_s=dt,
            mean_complete_step_wall_seconds=mean_cost,peak_process_working_set_bytes=m['peak_process_working_set_bytes']),
        geometry=dict(measurement_plane_z_m=p['layer'],**geometry,
            minimum_upstream_parent_state_center_gap_m=float(current_clearance[original_clearance>=0].min())),
        nominal_actuator_scale=dict(feed_m_s=p['feed'],time_to_advance_minimum_original_gap_s=nominal,
            meaning='NOMINAL_AXIAL_COMMAND_DISTANCE_ONLY_NOT_PREDICTED_MATERIAL_TRAVEL_OR_CROSSING_TIME'),
        fixed_step_linear_cost=dict(steps=steps,projected_wall_seconds=projection,
            projected_wall_days=projection/86400.,projected_wall_years=projection/(365.25*86400.),
            assumption='Current short-window average step cost held constant; NOT validated long-time forecast.'),
        decision=dict(repeat_short_steps_until_transport='NOT_AUTHORIZED_NOT_PROCESS_FEASIBLE',
            simply_relaxing_displacement_guard='INSUFFICIENT; NEED_PROCESS_DOMAIN_EDGE_CONTACT_AND_TIME_SCALE_VALIDATION',
            sustained_transport='NOT_ESTABLISHED',production_three_cases='NOT_RUN'))
    for rel,h in source_hashes.items():
        if sha(ROOT/rel)!=h:raise ValueError('input changed during read-only assessment')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh project-local report required')
    try:r=assess(args.case);code=0
    except Exception as error:r=dict(status='FAIL_FEASIBILITY_READBACK',error=repr(error),production_authorized=False);code=2
    out.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(r,allow_nan=False))
    return code


if __name__=='__main__':raise SystemExit(main())
