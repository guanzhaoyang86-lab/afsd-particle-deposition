"""Independent machine axes with exact phase endpoints, SI units throughout.

These are machine/actuator commands, never prescribed free-material paths.
Within each numerical interval every consumer uses the same linear-origin,
z-geodesic tool path. Phase boundaries must be actual integration boundaries.
"""
from dataclasses import dataclass
import copy
import math
import numpy as np

CLASSIFICATION = 'AFSD_INDEPENDENT_MACHINE_AXES_PHASE_PROGRAM_V1'
AXES = ('omega_rad_s', 'feed_m_s', 'traverse_m_s')


def finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('finite machine command required')
    return value


def rotation(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])


class ProcessProgram:
    def __init__(self, specification, shoulder_height_m):
        self.specification = copy.deepcopy(specification)
        self.layer = finite(shoulder_height_m)
        if self.layer <= 0 or not isinstance(specification, list) or not specification:
            raise ValueError('positive shoulder height and explicit phase list required')
        self.phases = []
        start = 0.; integral = np.zeros(3); previous = np.zeros(3); names = set()
        for entry in specification:
            name = entry['name']; duration = finite(entry['duration_s'])
            target = np.array([finite(entry[k]) for k in AXES])
            if (not isinstance(name, str) or not name or name in names or duration <= 0
                    or np.any(target < 0) or entry['transition'] not in ('smoothstep', 'hold')
                    or entry.get('status') != 'DECLARED_ENGINEERING_ASSUMPTION'
                    or not entry.get('source')):
                raise ValueError('explicit sourced nonnegative phase program required')
            if entry['transition'] == 'hold' and not np.array_equal(target, previous):
                raise ValueError('hold phase cannot introduce an instantaneous velocity jump')
            names.add(name)
            self.phases.append(dict(name=name, start=start, end=start+duration,
                duration=duration, initial=previous.copy(), target=target,
                integral0=integral.copy(), transition=entry['transition']))
            integral += duration*(previous+target)/2
            start += duration; previous = target
        self.end_s = start

    def _time(self, time_s):
        t = finite(time_s)
        if t < 0 or t > self.end_s:
            raise ValueError('time outside declared finite process program')
        return t

    def phase(self, time_s):
        t = self._time(time_s)
        return next((p for p in self.phases if t < p['end']), self.phases[-1])

    def integral(self, time_s):
        t = self._time(time_s); p = self.phase(t)
        local = t-p['start']; s = local/p['duration']
        return p['integral0']+p['initial']*local+(p['target']-p['initial'])*p['duration']*(s**3-.5*s**4)

    def velocity(self, time_s):
        t = self._time(time_s); p = self.phase(t)
        s = (t-p['start'])/p['duration']
        return p['initial']+(p['target']-p['initial'])*s*s*(3-2*s)

    def global_pose(self, time_s):
        angle, _, traverse = self.integral(time_s)
        return dict(origin=np.array([traverse, 0., self.layer]), rotation=rotation(angle))

    def bounded_dt(self, time_s, proposed_dt):
        t = self._time(time_s); dt = finite(proposed_dt)
        if dt <= 0 or t == self.end_s:
            raise ValueError('positive interval before program completion required')
        return min(dt, self.phase(t)['end']-t)

    def step(self, time_s, dt_s):
        return ProcessInterval(self, time_s, dt_s)

    def to_dict(self):
        return dict(classification=CLASSIFICATION, phases=copy.deepcopy(self.specification),
            shoulder_height_m=self.layer, end_s=self.end_s,
            free_material_trajectory_prescribed=False,
            within_step_path='LINEAR_WORLD_ORIGIN_Z_GEODESIC_LINEAR_AXIAL_ACTUATOR',
            startup_duration_source='DECLARED_ASSUMPTION_NOT_MEASURED_FROM_EDITED_VIDEO')


@dataclass(frozen=True, init=False)
class ProcessInterval:
    def __init__(self, history, t0, dt):
        if not isinstance(history, ProcessProgram):
            raise ValueError('explicit process program required')
        t0 = history._time(t0); dt = finite(dt); t1 = t0+dt
        if dt <= 0 or t1 <= t0 or t1 > history.phase(t0)['end']:
            raise ValueError('interval must end at or before the next phase boundary')
        start, end = history.integral(t0), history.integral(t1)
        if abs(end[0]-start[0]) >= math.pi:
            raise ValueError('interval rotation must be below pi radians')
        for k, v in dict(history=history, t0=t0, t1=t1, dt_s=dt,
                         start=start, end=end, phase_name=history.phase(t0)['name']).items():
            object.__setattr__(self, k, v)

    @property
    def feed_displacement_m(self):
        return -float(self.end[1]-self.start[1])

    def fraction(self, time_s):
        t = finite(time_s)
        if t < self.t0 or t > self.t1:
            raise ValueError('interval path cannot extrapolate')
        return (t-self.t0)/self.dt_s

    def pose(self, time_s):
        a = self.fraction(time_s)
        if time_s == self.t0:
            return self.history.global_pose(self.t0)
        if time_s == self.t1:
            return self.history.global_pose(self.t1)
        q = (1-a)*self.start+a*self.end
        return dict(origin=np.array([q[2], 0., self.history.layer]), rotation=rotation(q[0]))

    def prescribed_feed_at(self, time_s):
        return self.fraction(time_s)*self.feed_displacement_m

    def to_dict(self):
        a, b = self.pose(self.t0), self.pose(self.t1)
        return dict(classification=CLASSIFICATION, phase=self.phase_name,
            t0_s=self.t0, t1_s=self.t1, dt_s=self.dt_s,
            origin0_m=a['origin'].tolist(), origin1_m=b['origin'].tolist(),
            rotation0=a['rotation'].tolist(), rotation1=b['rotation'].tolist(),
            unwrapped_angle0_rad=float(self.start[0]), unwrapped_angle1_rad=float(self.end[0]),
            feed_displacement_m=self.feed_displacement_m,
            average_omega_rad_s=float((self.end[0]-self.start[0])/self.dt_s),
            average_traverse_m_s=float((self.end[2]-self.start[2])/self.dt_s),
            path='LINEAR_WORLD_ORIGIN_Z_GEODESIC_LINEAR_AXIAL_ACTUATOR')


def domain_and_inventory_audit(parameters, program):
    """Finite positive-x straight-pass coverage, not a convergence certificate."""
    p = parameters
    keys = ('substrate_width', 'substrate_depth', 'width', 'bar_length', 'dx', 'h_ratio',
            'tool_radius', 'bore_width', 'layer', 'actuator_depth', 'tool_top')
    if any(finite(p[k]) <= 0 for k in keys):
        raise ValueError('positive geometry required')
    if program.layer != p['layer']:
        raise ValueError('motion and geometry must use the same shoulder height')
    for key in ('substrate_width', 'substrate_depth', 'width', 'bar_length', 'actuator_depth'):
        cells = round(p[key]/p['dx'])
        if cells < 1 or abs(cells*p['dx']-p[key]) > 1e-12:
            raise ValueError('geometry must contain full particle cells: '+key)
    if p['width'] > p['bore_width'] or p['tool_radius'] <= p['bore_width']/math.sqrt(2):
        raise ValueError('open square bore and annular shoulder required')
    if p['tool_top'] <= p['layer'] or p['bar_length'] <= p['actuator_depth']+p['layer']:
        raise ValueError('finite tool and separated actuator zone required')
    angle, feed, distance = program.integral(program.end_s)
    buffer = 2*p['h_ratio']*p['dx']
    lateral = p['substrate_width']/2-p['tool_radius']
    forward = lateral-distance
    # Keep the entire material-labelled actuator cap above the outlet by half a cell.
    usable_feed = p['bar_length']-p['actuator_depth']-p['layer']-p['dx']/2
    failures = []
    if min(lateral, forward) < buffer-1e-14:
        failures.append('SUBSTRATE_DOES_NOT_COVER_FULL_PATH_PLUS_KERNEL_BUFFER')
    if feed > usable_feed+1e-14:
        failures.append('FINITE_FEEDSTOCK_ACTUATOR_RESERVE_EXHAUSTED')
    nsub = round(p['substrate_width']/p['dx'])**2*round(p['substrate_depth']/p['dx'])
    nbar = round(p['width']/p['dx'])**2*round(p['bar_length']/p['dx'])
    return dict(status='PASS_FINITE_DOMAIN_AND_INVENTORY_PLAN' if not failures else 'FAIL_DOMAIN_OR_INVENTORY_PLAN',
        passed=not failures, failures=failures, particle_count=nsub+nbar,
        substrate_count=nsub, bar_count=nbar, programmed_duration_s=program.end_s,
        programmed_feed_m=float(feed), programmed_traverse_m=float(distance),
        programmed_revolutions=float(angle/(2*math.pi)), required_edge_buffer_m=buffer,
        final_forward_edge_margin_m=float(forward), lateral_edge_margin_m=lateral,
        usable_axial_command_m=usable_feed, replenishment='FINITE_BAR_NO_AUTOMATIC_INSERTION',
        domain_convergence='NOT_TESTED', free_material_transport='NOT_PREDICTED_BY_COMMAND_INTEGRALS')
