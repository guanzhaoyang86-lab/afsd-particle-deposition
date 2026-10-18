"""Persistent identity lives in ParticleArray properties and follows all permutations."""
import numpy as np
from pysph.base.utils import get_particle_array

IDENTITY = ('particle_id', 'x0', 'y0', 'z0', 'material_id', 'component_id', 'tracer_type', 'wholebar_tag', 'actuator')


def make_array(x, component, p, first_id=1):
    x = np.asarray(x, float); n = len(x)
    pa = get_particle_array(name='solid', x=x[:, 0], y=x[:, 1], z=x[:, 2], h=np.full(n, p['dx']*p['h_ratio']), m=np.full(n, p['rho0']*p['dx']**3), rho=np.full(n, p['rho0']))
    for key in ('ax', 'ay', 'az', 'au', 'av', 'aw', 'arho', 'ae', 'e'):
        if key not in pa.properties: pa.add_property(key)
    edge = (component == 1) & (np.max(np.abs(x[:, :2]), axis=1) >= p['width']/2-p['dx']*.51)
    core = (component == 1) & (np.linalg.norm(x[:, :2], axis=1) <= p['dx']*.51)
    full = np.where(edge, 1, np.where(core, 2, 0))
    band = (x[:, 2] >= p['marker_low']) & (x[:, 2] <= p['marker_high'])
    actuator = (component == 1) & (x[:, 2] > p['bar_length']-p['actuator_depth'])
    props = dict(particle_id=np.arange(first_id, first_id+n, dtype=np.int64), material_id=np.ones(n, dtype=int), component_id=component, tracer_type=np.where(band, full, 0), wholebar_tag=full, actuator=actuator.astype(int))
    for key, values in props.items():
        pa.add_property(key, type='long', data=values)
    if not np.array_equal(pa.particle_id, props['particle_id']):
        raise OverflowError('persistent ID exceeds upstream integer storage')
    for axis, key in enumerate(('x0', 'y0', 'z0')): pa.add_property(key, data=x[:, axis])
    for key in ('ep', 'ep_rate', 'strain_rate', 'plastic_work', 'return_dissipation', 'feed_crossed'):
        pa.add_property(key)
    pa.add_property('temperature', data=np.full(n, p['temperature']))
    pa.add_property('S', stride=9, data=np.zeros((n, 9)))
    pa.add_property('D', stride=9, data=np.zeros((n, 9)))
    pa.u[actuator] = p['traverse']; pa.w[actuator] = -p['feed']
    pa.set_output_arrays(list(pa.properties))
    return pa


def create_geometry(p):
    dx = p['dx']
    def centers(width):
        count = int(round(width/dx))
        if abs(count*dx-width) > 1e-12: raise ValueError('geometry must fit full particle cells')
        return (np.arange(count)+.5)*dx-width/2
    xy = centers(p['substrate_width'])
    z = centers(p['substrate_depth'])-p['substrate_depth']/2
    substrate = np.array(np.meshgrid(xy, xy, z, indexing='ij')).reshape(3, -1).T
    barxy = centers(p['width']); barz = centers(p['bar_length'])+p['bar_length']/2
    bar = np.array(np.meshgrid(barxy, barxy, barz, indexing='ij')).reshape(3, -1).T
    return make_array(np.vstack((substrate, bar)), np.r_[np.zeros(len(substrate), int), np.ones(len(bar), int)], p)


def identity_map(pa):
    return {int(pa.particle_id[i]): tuple(getattr(pa, k)[i].item() for k in IDENTITY[1:]) for i in range(pa.get_number_of_particles())}


def append_feedstock(pa, x, p, next_id):
    """Tested identity interface, not enabled as an infinite inlet by the runner."""
    if next_id <= int(pa.particle_id.max()): raise ValueError('ID reuse forbidden')
    incoming = make_array(x, np.ones(len(x), int), p, next_id)
    pa.append_parray(incoming)
    return next_id+len(x)
