"""SI finite flat-shoulder/open rotating-square geometry, no particle animation."""
import numpy as np


def tool_pose(t, traverse, omega, layer):
    return np.array([traverse*t, 0., layer]), omega*t


def local_xy(x, t, traverse, omega):
    angle = omega*t
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    return (x[:, :2]-[traverse*t, 0])@rotation, rotation


def as_rs(y, omega, traverse):
    """+z CCW, traverse +x: local rotation v_x=-omega*y; AS has same sign."""
    product = -omega*np.asarray(y)*traverse
    return np.where(product > 0, 1, np.where(product < 0, -1, 0))


def assert_outside_tool(x, t, p):
    """Fail closed when finite-face broad-phase would miss a deeply buried point.

    The solid is the finite circular prism minus the rotating square bore.
    A particle centre entering this solid already exceeds the supported 0.1dx
    overlap gate (particle radius=0.5dx). This check never projects particles.
    Top/outer-edge contact remains unsupported; this guard is not its substitute.
    """
    xy, _ = local_xy(x, t, p['traverse'], p['omega'])
    inside = ((x[:, 2] > p['layer']) & (x[:, 2] < p['tool_top']) &
              (np.linalg.norm(xy, axis=1) < p['tool_radius']) &
              (np.max(np.abs(xy), axis=1) > p['bore_width']/2))
    if np.any(inside):
        raise ValueError('particle centre entered rigid tool solid; finite-face contact lost: indices='+
                         repr(np.flatnonzero(inside)[:20].tolist()))


def surfaces(x, t, p):
    """Each row: indices, signed gap, inward normal, wall velocity, wall T.

    Finite analytic faces are a DEBUG alternative to the author's triangle mesh.
    Contact radius is dx/2; square corners may have two unilateral constraints.
    """
    rad = p['dx']/2
    xy, rotation = local_xy(x, t, p['traverse'], p['omega'])
    wall_v = np.column_stack((p['traverse']-p['omega']*x[:, 1], p['omega']*(x[:, 0]-p['traverse']*t), np.zeros(len(x))))
    # Rigid base.
    ids = np.flatnonzero(x[:, 2]+p['substrate_depth'] < rad)
    yield ids, x[ids, 2]+p['substrate_depth'], np.tile([0., 0., 1.], (len(ids), 1)), np.zeros((len(ids), 3)), p['base_temperature']
    # Flat shoulder annulus with open square bore (no pin / solid cylinder).
    annulus = (np.max(np.abs(xy), axis=1) >= p['bore_width']/2) & (np.linalg.norm(xy, axis=1) <= p['tool_radius'])
    ids = np.flatnonzero(annulus & (x[:, 2] < p['layer']+rad) & (p['layer']-x[:, 2] < rad))
    yield ids, p['layer']-x[ids, 2], np.tile([0., 0., -1.], (len(ids), 1)), wall_v[ids], p['tool_temperature']
    for axis in (0, 1):
        for sign in (-1, 1):
            distance = p['bore_width']/2-sign*xy[:, axis]
            mask = (x[:, 2] >= p['layer']) & (x[:, 2] <= p['tool_top']) & (np.abs(xy[:, 1-axis]) < p['bore_width']/2+rad) & (distance < rad) & (distance > -2*rad)
            ids = np.flatnonzero(mask)
            local_n = np.zeros(2); local_n[axis] = -sign
            normal = np.r_[rotation@local_n, 0.]
            yield ids, distance[ids], np.tile(normal, (len(ids), 1)), wall_v[ids], p['tool_temperature']
