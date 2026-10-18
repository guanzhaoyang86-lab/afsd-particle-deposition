"""Explicit hot-start boundary history, independent of every old process state.

Nominal omega/traverse/feed are multiplied by f(s)=3*s**2-2*s**3 for
s=t/T in [0,1], then by one. Its exact analytical integral is used for
every interval endpoint. Within ONE step the declared numerical path has
linear world translation, a constant-rate z geodesic rotation and linear
prescribed actuator translation joining these endpoints. It is intentionally
not the curved global ramp at intermediate times. The interval object is the
single source for contact, support/shell certificates and passive observers.

This is a new startup boundary-condition assumption, not a temporal/mass/
modulus scaling, an inertial restart or an assertion of improved convergence.
"""
from dataclasses import dataclass
import math
import numpy as np

RAMP_DURATION_S=.020
MOTION_CLASSIFICATION='SMOOTHSTEP_VELOCITY_EXACT_ENDPOINTS_PIECEWISE_LINEAR_GEODESIC_STEP_V1'


def _finite(value,name):
    value=float(value)
    if not math.isfinite(value):raise ValueError('finite '+name+' required')
    return value


def _rotation(angle):
    c,s=math.cos(angle),math.sin(angle)
    return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])


@dataclass(frozen=True,init=False)
class RampedHotMotion:
    omega: float
    traverse: float
    feed: float
    layer: float
    ramp_duration_s: float

    def __init__(self,parameters,ramp_duration_s=RAMP_DURATION_S):
        for name in ('omega','traverse','feed','layer'):
            object.__setattr__(self,name,_finite(parameters[name],name))
        duration=_finite(ramp_duration_s,'ramp duration')
        if duration<=0:raise ValueError('positive ramp duration required')
        object.__setattr__(self,'ramp_duration_s',duration)

    def _time(self,t):
        t=_finite(t,'nonnegative motion time')
        if t<0:raise ValueError('independent motion time cannot precede zero')
        return t

    def velocity_factor(self,t):
        t=self._time(t)
        if t>=self.ramp_duration_s:return 1.
        s=t/self.ramp_duration_s
        return s*s*(3.-2.*s)

    def exposure(self,t):
        """Analytical integral of the velocity factor, in seconds."""
        t=self._time(t);T=self.ramp_duration_s
        if t>=T:return t-.5*T
        s=t/T
        return T*s*s*s*(1.-.5*s)

    def exposure_increment(self,t0,t1):
        """Factored integral difference avoids subtracting nearby quartics."""
        a,b=self._time(t0),self._time(t1)
        if b<a:raise ValueError('motion interval is reversed')
        if b==a:return 0.
        T=self.ramp_duration_s
        if a>=T:return b-a
        if b>T:return self.exposure_increment(a,T)+(b-T)
        x,y=a/T,b/T
        return (b-a)*((y*y+x*y+x*x)-.5*(y*y*y+x*y*y+x*x*y+x*x*x))

    def global_pose(self,t):
        """Exact global ramp endpoint; not the declared path inside a step."""
        e=self.exposure(t)
        return dict(origin=np.array([self.traverse*e,0.,self.layer]),rotation=_rotation(self.omega*e))

    def to_dict(self):
        return dict(classification=MOTION_CLASSIFICATION,ramp_duration_s=self.ramp_duration_s,
            velocity_factor='3*s^2-2*s^3 for 0<=s=t/T<=1; 1 thereafter',
            integrated_factor='T*(s^3-s^4/2) until T; t-T/2 thereafter',
            nominal_omega_rad_s=self.omega,nominal_traverse_m_s=self.traverse,
            nominal_feed_m_s=self.feed,shoulder_height_m=self.layer,
            step_path='LINEAR_WORLD_ORIGIN_Z_GEODESIC_ROTATION_LINEAR_ACTUATOR',
            interpretation='NEW_STARTUP_BOUNDARY_HISTORY_ASSUMPTION_NOT_PHYSICAL_VALIDATION')

    def step(self,t0,dt):return IntervalMotion(self,t0,dt)


@dataclass(frozen=True,init=False)
class IntervalMotion:
    history: RampedHotMotion
    t0: float
    t1: float
    dt_s: float
    exposure0_s: float
    exposure1_s: float
    exposure_delta_s: float

    def __init__(self,history,t0,dt):
        if not isinstance(history,RampedHotMotion):raise ValueError('explicit ramp history required')
        t0=history._time(t0);dt=_finite(dt,'positive interval')
        if dt<=0 or not math.isfinite(t0+dt) or t0+dt<=t0:raise ValueError('resolvable positive motion interval required')
        t1=t0+dt;e0=history.exposure(t0);e1=history.exposure(t1)
        if abs(history.omega*(e1-e0))>=math.pi:
            raise ValueError('single step must have an unambiguous z geodesic below pi radians')
        for key,value in dict(history=history,t0=t0,t1=t1,dt_s=dt,exposure0_s=e0,
                exposure1_s=e1,exposure_delta_s=history.exposure_increment(t0,t1)).items():
            object.__setattr__(self,key,value)

    @property
    def feed_displacement_m(self):return -self.history.feed*self.exposure_delta_s

    def fraction(self,t):
        t=_finite(t,'interval motion time')
        if t<self.t0 or t>self.t1:raise ValueError('interval pose cannot extrapolate outside its declared endpoints')
        if t==self.t0:return 0.
        if t==self.t1:return 1.
        return (t-self.t0)/(self.t1-self.t0)

    def pose(self,t):
        alpha=self.fraction(t)
        # Return exact global endpoints, and interpolate only the declared step path.
        if alpha==0.:return self.history.global_pose(self.t0)
        if alpha==1.:return self.history.global_pose(self.t1)
        exposure=(1.-alpha)*self.exposure0_s+alpha*self.exposure1_s
        return dict(origin=np.array([self.history.traverse*exposure,0.,self.history.layer]),
            rotation=_rotation(self.history.omega*exposure))

    def prescribed_feed_at(self,t):
        """Actuator displacement relative to this interval's start."""
        return self.fraction(t)*self.feed_displacement_m

    def to_dict(self):
        a,b=self.pose(self.t0),self.pose(self.t1)
        return dict(classification=MOTION_CLASSIFICATION,t0_s=self.t0,t1_s=self.t1,dt_s=self.dt_s,
            exposure0_s=self.exposure0_s,exposure1_s=self.exposure1_s,exposure_increment_s=self.exposure_delta_s,
            origin0_m=a['origin'].tolist(),origin1_m=b['origin'].tolist(),
            rotation0=a['rotation'].tolist(),rotation1=b['rotation'].tolist(),
            unwrapped_angle0_rad=self.history.omega*self.exposure0_s,
            unwrapped_angle1_rad=self.history.omega*self.exposure1_s,
            feed_displacement_m=self.feed_displacement_m,
            nominal_velocity_factor0=self.history.velocity_factor(self.t0),
            nominal_velocity_factor1=self.history.velocity_factor(self.t1),
            actual_step_average_velocity_factor=self.exposure_delta_s/self.dt_s,
            path='LINEAR_WORLD_ORIGIN_Z_GEODESIC_ROTATION_LINEAR_ACTUATOR',
            global_ramp_not_used_inside_step=True)
