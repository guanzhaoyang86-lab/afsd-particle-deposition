"""Initial-material axial actuator, never a free-material trajectory driver."""
import numpy as np


def control(z0,time_s,base_z,height,depth,speed,release_time):
    z0=np.asarray(z0,float)
    if z0.ndim!=1 or not np.isfinite(z0).all():raise ValueError('finite original z vector required')
    if not np.isfinite([time_s,base_z,height,depth,speed,release_time]).all():raise ValueError('finite actuator inputs required')
    if time_s<0 or height<=0 or not 0<depth<height or speed<=0 or release_time<=0:raise ValueError('invalid actuator program')
    eligible=(z0>base_z+height-depth)&(z0<=base_z+height)
    active=eligible&(time_s<release_time)
    free=np.ones((len(z0),3));free[active,2]=0.
    prescribed=np.zeros((len(z0),3));prescribed[active,2]=-speed
    return dict(eligible=eligible,active=active,free_axes=free,prescribed_velocity=prescribed,
        phase='DRIVEN' if time_s<release_time else 'RELEASED',
        commanded_displacement_m=-speed*min(time_s,release_time),
        displacement_status='ACTUATOR_COMMAND_ONLY_NOT_FREE_MATERIAL_PATH')


def validate_interval(time_s,dt,release_time):
    if not np.isfinite([time_s,dt,release_time]).all() or time_s<0 or dt<=0 or release_time<=0:
        raise ValueError('invalid actuator interval')
    # The companion only registers dt values aligning the discontinuity. This
    # guard forbids silently straddling the release within a frozen-constraint kick.
    if time_s+1e-20<release_time<time_s+dt-1e-20:
        raise ValueError('actuator release must lie on a step boundary')
